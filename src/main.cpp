
#include <iostream>

#ifdef  _WIN32
#include "win\dirent.h"
#else
#include "dirent.h"
#endif

#include "WfipsGrid.h"
#include "ogrsf_frmts.h"
#include "gdal_alg.h"

#include <omp.h>
#include <set>
#include <vector>
#include <unordered_map>
#include <cstdio>
#include <gdal_utils.h>

#include "sqlite3.h"

#include "WfipsGridData.h"

using std::multimap;
using std::string;
using std::vector;
//using std::printf;

/***************************************************************************
// Data Structures
***************************************************************************/

//data structure to hold bounding box
struct SBoundingBox
{
    double maxX;
    double maxY;
    double minX;
    double minY;
};

struct WfipsData
{
    OGRSpatialReference spatialReference;
    //vector<OGRPoint> cellCentroids;
    //vector<SBoundingBox> cellBoundingBoxes;
    CWfipsGridData gridData;

    double GeoTransform[6];

    int numRows = -1;
    int numCols = -1;
};

struct FireshedData
{
    //vector<int> fireNumbers;
    vector<vector<int>> fireOriginCells;
    vector<SBoundingBox> cellBoundingBoxes;
    vector<vector<unordered_map<int, int>>> wfipscellsToFireOriginsForSingleFile;
    vector<vector<int>> originCellsForWfipscell;
    unordered_map<int, int> finalIndexToWfipsCell;
    vector<vector<int>> frequenciesForWfipscell;
};

bool BoundingBoxCheck(SBoundingBox innerBoundingBox, struct SBoundingBox outerBoundingBox);

void ReadFromShapeFileToMemory(bool verbose, const int shapeFileListSize, const int shapeFileIndex, const std::string shapeFilePath, const std::string shapeFileName, FireshedData& fireshedData, WfipsData& wfipsData);
void FillOriginDataForSingleFire(vector<int>& rasterBuffer, const int shapeFileIndex, const int fireIndex, const OGREnvelope fireBoundingBox, FireshedData& fireshedData, WfipsData& wfipsData);
void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData);
int FillWfipsData(WfipsData& wfipsData, std::string dataPath);
void CreateFireShedDB(string outPath, FireshedData& fireshedData, WfipsData& wfipsData);

static const double cellHalfWidth = 1000; // 1 km

int main(int argc, char *argv[])
{
    DIR *dir;
    struct dirent *ent;

    string fileName = "";
    string extension = "";
    size_t pos = -1;
    int type = 0;
    string dataPath = "";
    string outPath = "";
    int rc = 0;
    bool verbose = false;

    if (argc == 3)
    {
        dataPath = argv[1];
        outPath = argv[2];
    }
    else if (argc == 4)
    {
        dataPath = argv[1];
        outPath = argv[2];
        string verboseTest = argv[3];
        if (verboseTest == "verbose")
        {
            verbose = true;
        }
        else
        {
            printf("Error: need path to shapefiles and db output directory as arguments");
            printf("\n    for progress info on console, enter \"verbose\" as third argument");
            return EXIT_FAILURE;
        }
    }
    else
    {
        printf("Error: need path to shapefiles and db output directory as arguments");
        printf("\nfor progress info on console, enter \"verbose\" as third argument");
        return EXIT_FAILURE;
    }

    if ((dataPath.back() != '/') && (dataPath.back() != '\\'))
    {
#ifdef WIN32
        dataPath.push_back('\\');
        outPath.push_back('\\');
#else
        dataPath.push_back('/');
        outPath.push_back('/');
#endif
    }

    vector<string> shapefileNameList;
    vector<string> shapefilePathList;

    if ((dir = opendir(argv[1])) != nullptr)
    {
        /* readt all the files and directories within directory */
        while ((ent = readdir(dir)) != nullptr)
        {
            fileName = ent->d_name;

            pos = fileName.find_last_of('.');
            if (pos != string::npos)
            {
                extension = fileName.substr(pos, fileName.size());
                if (extension == ".SHP" || extension == ".shp")
                {
                    //printf("%s\n", temp.c_str());
                    shapefileNameList.push_back(fileName);
                    shapefilePathList.push_back(dataPath + fileName);
                }
            }
        }
        closedir(dir);
    }
    else
    {
        /* could not open directory */
        perror("");
        return EXIT_FAILURE;
    }

    if (shapefilePathList.size() < 1)
    {
        printf("Error: No shapefiles found, exiting program\n");
        return EXIT_FAILURE;
    }

    const int shapefileListSize = shapefilePathList.size();

    GDALAllRegister();

    char *MapESRIProjStrings[] =
    {
        "",
        "PROJCS[\"Albers\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Albers\"],PARAMETER[\"standard_parallel_1\",29.5],PARAMETER[\"standard_parallel_2\",45.5],PARAMETER[\"latitude_of_origin\",23],PARAMETER[\"central_meridian\",-96],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"METERS\",1]]",
        "PROJCS[\"WGS 84 / Pseudo - Mercator\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Mercator\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]"
    };

    CPLStringList papszPrj5070;
    papszPrj5070.AddString(MapESRIProjStrings[1]);

    WfipsData wfipsData;
    rc = wfipsData.spatialReference.importFromESRI(papszPrj5070);
    rc = FillWfipsData(wfipsData, dataPath);
    if (rc != 0)
    {
        printf("Error: WFIPS data loading failed, exiting program\n");
        return EXIT_FAILURE;
    }
    printf("WFIPS data populated\n");

    FireshedData fireshedData;

    fireshedData.wfipscellsToFireOriginsForSingleFile.resize(shapefileListSize);
    fireshedData.fireOriginCells.resize(shapefileListSize);
    fireshedData.wfipscellsToFireOriginsForSingleFile.resize(shapefileListSize);

    double currentProgress = 0;
    for (int shapefileIndex = 0; shapefileIndex < shapefileListSize; shapefileIndex++)
    {
        ReadFromShapeFileToMemory(verbose, shapefileListSize, shapefileIndex, dataPath, shapefileNameList[shapefileIndex], fireshedData, wfipsData);
        currentProgress = (shapefileIndex / (shapefileListSize * 1.0)) * 100.00;
        if (verbose)
        {
            printf("processed %d files out of %d in\n", shapefileIndex, shapefileListSize);
            printf("    %s\n    %4.2f percent of all files to be processed are complete\n\n", dataPath.c_str(), currentProgress);
        }
    }

    printf("processed %d files out of %d in\n", shapefileListSize, shapefileListSize);
    printf("    %s\n    100 percent of all files to be processed are complete\n\n", dataPath.c_str());

    ConsolidateFinalData(shapefileListSize, fireshedData);
    CreateFireShedDB(outPath, fireshedData, wfipsData);

    return 0;
}

/***************************************************************************
// Shapefile Reading Function
***************************************************************************/
void ReadFromShapeFileToMemory(bool verbose, const int shapefileListSize, const int shapefileIndex, const std::string shapefilePath, const std::string shapefileName, FireshedData& fireshedData, WfipsData& wfipsData)
{
    OGRErr error;
    string shapefileFullPath = shapefilePath + shapefileName;
    GDALDataset *poShapefileDS = static_cast<GDALDataset*>(GDALOpenEx(shapefileFullPath.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL));

    char *MapESRIProjStrings[] =
    {
        "",
        "PROJCS[\"Albers\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Albers\"],PARAMETER[\"standard_parallel_1\",29.5],PARAMETER[\"standard_parallel_2\",45.5],PARAMETER[\"latitude_of_origin\",23],PARAMETER[\"central_meridian\",-96],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"METERS\",1]]",
        "PROJCS[\"WGS 84 / Pseudo - Mercator\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Mercator\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]"
    };

    OGRLayer  *poLayer = poShapefileDS->GetLayer(0);
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();

    OGRwkbGeometryType LayerGeometryType = poLayer->GetGeomType();
    int NumberOfFeatures = poLayer->GetFeatureCount(true);

    poLayer->ResetReading();

    int sizeInAcres = 0;
    int fireNumber;

    enum
    {
        fire_number = 0,
        acres = 7,
        x_val = 8,
        y_val = 9
    };

    double currentProgress = 0;

    int nBufXSize = wfipsData.gridData.GetNumX();
    int nBufYSize = wfipsData.gridData.GetNumY();
    int nBandCount = 1;
    GDALDataType eType = GDT_Int32;
    int nDataTypeSize = GDALGetDataTypeSizeBytes(eType);
    int bandList[1] = { 1 };
    std::vector<double> geomBurnValue(nBandCount, 1.0);
    std::vector<OGRGeometryH> ogrBurnGeometries;

    if (wkbFlatten(LayerGeometryType) == wkbPolygon)
    {
        //Read Polygon Shapefile into vector of Geometries
        for (int geometryIndex = 0; geometryIndex < NumberOfFeatures; geometryIndex++)
        {
            OGRFeature *poInputFeature;
            poInputFeature = poLayer->GetNextFeature();
            OGRGeometry *poGeometry;

            poGeometry = poInputFeature->GetGeometryRef();

            double x = 0.0;
            double y = 0.0;

            for (int iField = 0; iField < 10; iField++)
            {
                OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);

                if (iField == fire_number)
                {
                    fireNumber = poInputFeature->GetFieldAsInteger(iField);
                }
                if (iField == acres)
                {
                    sizeInAcres = poInputFeature->GetFieldAsInteger(iField);
                }
                else if (iField == x_val)
                {
                    x = poInputFeature->GetFieldAsDouble(iField);
                }
                else if (iField == y_val)
                {
                    y = poInputFeature->GetFieldAsDouble(iField);
                }
            }

            if (sizeInAcres >= 150)
            {
                int originCell = (wfipsData.gridData.WG_GetCellIndex(x, y));
                fireshedData.fireOriginCells[shapefileIndex].push_back(originCell);
                poGeometry->closeRings();
                ogrBurnGeometries.push_back((OGRGeometryH)poGeometry->clone());
            }
            OGRFeature::DestroyFeature(poInputFeature);
        }
    }
    GDALClose(poShapefileDS);

    int firesProcessed = 0;
    const int totalFires = ogrBurnGeometries.size();

    fireshedData.wfipscellsToFireOriginsForSingleFile[shapefileIndex].resize(totalFires);
    fireshedData.fireOriginCells[shapefileIndex].resize(totalFires);

    int numCores = omp_get_num_procs();
    omp_set_num_threads(numCores);

    vector<vector<int>> buffer;
    vector<GDALDataset*> inMemoryRaster;
    GDALDriverH hMemDriver = GDALGetDriverByName("MEM");

    buffer.resize(numCores);
    inMemoryRaster.resize(numCores);

    // Parallel
#pragma omp parallel for shared(inMemoryRaster)
    for (int threadIndex = 0; threadIndex < numCores; threadIndex++)
    {
        // Create in-memory rasters
        GDALDatasetH hRasterDS = GDALCreate(hMemDriver, "", nBufXSize, nBufYSize, nBandCount, eType, NULL);
        inMemoryRaster[threadIndex] = (GDALDataset *)hRasterDS;
        buffer[threadIndex].resize(nBufXSize*nBufYSize);
    }

#pragma omp parallel for shared(firesProcessed, currentProgress)
    for (int fireIndex = 0; fireIndex < totalFires; fireIndex++)
    {
        int threadID = omp_get_thread_num();

        // Assign GeoTransform parameters
        GDALDatasetH hRasterDS = inMemoryRaster[threadID];
        OGRErr ogrError = GDALSetGeoTransform(hRasterDS, wfipsData.GeoTransform);
        ogrError = GDALSetProjection(hRasterDS, MapESRIProjStrings[1]);

        OGRGeometry *poGeometry = (OGRGeometry *)ogrBurnGeometries[fireIndex];

        OGREnvelope envelope;
        poGeometry->getEnvelope(&envelope);

        char** options = nullptr;
        options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");

        CPLErr err = CE_None;

        err = GDALRasterizeGeometries(hRasterDS, nBandCount, bandList,
            1, &(ogrBurnGeometries[fireIndex]), NULL, NULL,
            geomBurnValue.data(), options, NULL, NULL);

        //// Test
        //if (err == CE_None)
        //{
        //    //Create Tiff from in-memory raster
        //    GDALDatasetH hDstDS;
        //    string geoTiffDriverName = "GTiff";

        //    GDALDriverH hGeoTiffDriver = GDALGetDriverByName(geoTiffDriverName.c_str());
        //    if (hGeoTiffDriver == NULL)
        //    {
        //        printf("%s driver not available.\n", geoTiffDriverName.c_str());
        //        exit(1);
        //    }
        //    string destination = shapeFilePath + "\\test\\test.tif";
        //    hDstDS = GDALCreateCopy(hGeoTiffDriver, destination.c_str(), hMemDset, FALSE,
        //        NULL, NULL, NULL);
        //    /* Once we're done, close properly the dataset */
        //    if (hDstDS != NULL)
        //    {
        //        GDALClose(hDstDS);
        //    }
        //}
        GDALDataset* rasterDS = inMemoryRaster[threadID];

        rasterDS->RasterIO(GF_Read, 0, 0, nBufXSize, nBufYSize, buffer[threadID].data(), nBufXSize, nBufYSize, GDT_Int32, 1, NULL, 0, 0, NULL);
        FillOriginDataForSingleFire(buffer[threadID], shapefileIndex, fireIndex, envelope, fireshedData, wfipsData);

        std::fill(buffer[threadID].data(), buffer[threadID].data(), 0);

        rasterDS->RasterIO(GF_Write, 0, 0, nBufXSize, nBufYSize, buffer[threadID].data(), nBufXSize, nBufYSize, GDT_Int32, 1, NULL, 0, 0, NULL);

#pragma omp atomic
        firesProcessed++;

        if (verbose && (firesProcessed % 1000) == 0)
        {
            //Critical Section
#pragma omp critical
            {
                currentProgress = (firesProcessed / (totalFires * 1.0)) * 100.00;
                printf("processed %d fires out of %d in file\n", firesProcessed, totalFires);
                printf("    %s\n    %4.2f percent of current file is complete\n\n", shapefileName.c_str(), currentProgress);
                printf("processed %d files out of %d in\n", shapefileIndex, shapefileListSize);
                currentProgress = (shapefileIndex / (shapefileListSize * 1.0)) * 100.00;
                printf("    %s\n    %4.2f percent of all files to be processed are complete\n\n", shapefilePath.c_str(), currentProgress);
            }
            //End Critical Section
        }
    }

#pragma omp parallel for shared(inMemoryRaster)
    for (int threadIndex = 0; threadIndex < numCores; threadIndex++)
    {
        // Destroy in-memory rasters
        GDALClose(inMemoryRaster[threadIndex]);
        inMemoryRaster[threadIndex] = NULL;
    }
    // End Parallel

    inMemoryRaster.clear();

    printf("file\n    %s\n   is 100 percent complete\n", shapefileName.c_str());
}

bool BoundingBoxCheck(SBoundingBox innerBoundingBox, struct SBoundingBox outerBoundingBox)
{
    return ((innerBoundingBox.minX >= outerBoundingBox.minX) && (innerBoundingBox.maxX <= outerBoundingBox.maxX) && (innerBoundingBox.minY >= outerBoundingBox.minY) && (innerBoundingBox.maxY <= outerBoundingBox.maxY));
}

int FillWfipsData(WfipsData& wfipsData, std::string dataPath)
{
    int rc = wfipsData.gridData.LoadData((const char*)dataPath.c_str());
    std::printf("WFIPS grid data initialized\n");

    wfipsData.numRows = wfipsData.gridData.GetNumY();
    wfipsData.numCols = wfipsData.gridData.GetNumX();

    wfipsData.gridData.GetWfipsGridGeotransform(wfipsData.GeoTransform);

    return rc;
}

void FillOriginDataForSingleFire(vector<int>& rasterBuffer, const int shapeFileIndex, const int geometryIndex, const OGREnvelope fireBoundingBox, FireshedData& fireshedData, WfipsData& wfipsData)
{
    int origin = fireshedData.fireOriginCells[shapeFileIndex][geometryIndex];

    double maxX = fireBoundingBox.MaxX;
    double maxY = fireBoundingBox.MaxY;
    double minX = fireBoundingBox.MinX;
    double minY = fireBoundingBox.MinY;

    int nBufXSize = wfipsData.gridData.GetNumX();
    int nBufYSize = wfipsData.gridData.GetNumY();

    int upperLeft = wfipsData.gridData.WG_GetCellIndex(minX, maxY);
    int upperRight = wfipsData.gridData.WG_GetCellIndex(maxX, maxY);
    int lowerLeft = wfipsData.gridData.WG_GetCellIndex(minX, minY);
    int lowerRight = wfipsData.gridData.WG_GetCellIndex(maxX, minY);

    int wfipsCell;
    int minRow;
    int maxRow;
    int minCol;
    int maxCol;

    wfipsData.gridData.GetRowColFromWfipsCell(upperLeft, &minRow, &minCol);
    wfipsData.gridData.GetRowColFromWfipsCell(lowerRight, &maxRow, &maxCol);

    // Get window for reading raster buffer, add a 1 cell buffer around perimeter if possible
    if (maxRow < nBufXSize)
    {
        maxRow++;
    }
    if (minRow > 0)
    {
        minRow--;
    }

    if (maxCol < nBufYSize)
    {
        maxCol++;
    }
    if (minCol > 0)
    {
        minCol--;
    }

    for (int row = minRow; row < maxRow; row++)
    {
        for (int col = minCol; col < maxCol; col++)
        {
            wfipsCell = wfipsData.gridData.GetWfipsCellIndex(row, col);
            if (rasterBuffer[wfipsCell] == 1)
            {
                fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex][geometryIndex].insert(std::make_pair(wfipsCell, origin));
            }
        }
    }
}

void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData)
{
    multiset<pair<int, int>> totalWfipscellsToFireOrigins;

    int wfipscell = -1;
    int origin = -1;
    int wfipscellPrevious = -1;
    int originPrevious = -1;
    int wfipsCellIndex = -1;
    int frequency = 1;

    for (int shapeFileIndex = 0; shapeFileIndex < num_shape_files; shapeFileIndex++)
    {
        int numFires = fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].size();
        for (int fireIndex = 0; fireIndex < numFires; fireIndex++)
        {
            for (auto iterator = fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex][fireIndex].begin(); iterator != fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex][fireIndex].end(); iterator++)
            {
                wfipscell = iterator->first;
                origin = iterator->second;
                totalWfipscellsToFireOrigins.insert(std::make_pair(wfipscell, origin));
            }
        }
        fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].clear();
    }

    for (auto iterator = totalWfipscellsToFireOrigins.begin(); iterator != totalWfipscellsToFireOrigins.end(); iterator++)
    {
        wfipscell = iterator->first;
        origin = iterator->second;
        bool wfipsCellChanged = false;
        bool originChanged = false;

        if (wfipscell != wfipscellPrevious) // Wfispcell changed, add new vector rows
        {
            wfipsCellChanged = true;
            wfipscellPrevious = wfipscell;
            wfipsCellIndex++;
            // Add current wfipscell index to map so it can be retrieved by row index later
            fireshedData.finalIndexToWfipsCell.insert(std::make_pair(wfipsCellIndex, wfipscell));
            vector<int> originVectorRow;
            fireshedData.originCellsForWfipscell.push_back(originVectorRow);
            vector<int> frequencyVectorRow;
            fireshedData.frequenciesForWfipscell.push_back(frequencyVectorRow);
        }
        if (originPrevious != origin) // Origin changed
        {
            originChanged = true;
            originPrevious = origin;
        }

        if (originChanged || wfipsCellChanged) // Add element to current vector row
        {
            frequency = 1;
            fireshedData.originCellsForWfipscell[wfipsCellIndex].push_back(origin);
            fireshedData.frequenciesForWfipscell[wfipsCellIndex].push_back(frequency);
        }
        else // Increment frequency and overwrite element containing frequency count for current wfipscell
        {
            frequency++;
            fireshedData.frequenciesForWfipscell[wfipsCellIndex].back() = frequency;
        }
    }

    totalWfipscellsToFireOrigins.clear();
}

void CreateFireShedDB(string outPath, FireshedData& fireshedData, WfipsData& wfipsData)
{
    int rc = 0;
    sqlite3 *db = NULL;
    sqlite3_stmt *stmt;

    char *sqlErrMsg = 0;

    string createFireshedDBSQLString = "CREATE TABLE IF NOT EXISTS firesheds(" \
        "wfipscell INTEGER," \
        "origin INTEGER,"
        "frequency INTEGER, " \
        "x REAL," \
        "y REAL)";

    string insertSQLString = "INSERT INTO firesheds(" \
        "wfipscell, " \
        "origin, " \
        "frequency, " \
        "x, " \
        "y) " \
        "VALUES(" \
        ":wfipscell, " \
        ":origin, " \
        ":frequency, " \
        ":x, " \
        ":y)";

    struct ColumnIndex
    {
        enum ColumnIndexEnum
        {
            wfipscell = 1,
            origin = 2,
            frequency = 3,
            x = 4,
            y = 5
        };
    };

    int wfipscell = -1;
    int bindColumnIndex = -1;
    int frequency = -1;
    int origin = -1;
    int err = -1;
    double cellCenterX = -1;
    double cellCenterY = -1;
    double cellMinX = -1;
    double cellMinY = -1;
    double cellMaxX = -1;
    double cellMaxY = -1;

    string outPutDbFullPath = outPath + "firesheds.db";
    rc = sqlite3_open_v2(outPutDbFullPath.c_str(), &db,
        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
        NULL);

    rc = sqlite3_exec(db, "DROP TABLE IF EXISTS 'firesheds'", NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, createFireshedDBSQLString.c_str(), NULL, NULL, &sqlErrMsg);

    rc = sqlite3_prepare_v2(db, insertSQLString.c_str(), -1, &stmt, NULL); // Prepare SQL statement

    rc = sqlite3_exec(db, "PRAGMA SYNCHRONOUS=OFF", NULL, NULL, NULL);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &sqlErrMsg);

    for (int wfipsCellIndex = 0; wfipsCellIndex < fireshedData.originCellsForWfipscell.size(); wfipsCellIndex++)
    {
        for (int originCellIndex = 0; originCellIndex < fireshedData.originCellsForWfipscell[wfipsCellIndex].size(); originCellIndex++)
        {
            wfipscell = fireshedData.finalIndexToWfipsCell.at(wfipsCellIndex);
            bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":wfipscell");
            rc = sqlite3_bind_int(stmt, bindColumnIndex, wfipscell);

            origin = fireshedData.originCellsForWfipscell[wfipsCellIndex][originCellIndex];
            bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":origin");
            rc = sqlite3_bind_int(stmt, bindColumnIndex, origin);

            frequency = fireshedData.frequenciesForWfipscell[wfipsCellIndex][originCellIndex];
            bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":frequency");
            rc = sqlite3_bind_int(stmt, bindColumnIndex, frequency);

            wfipsData.gridData.WG_GetCellCoords(wfipscell, &cellMinX, &cellMinY, &cellMaxX, &cellMaxY);

            cellCenterX = cellMinX + cellHalfWidth;
            cellCenterY = cellMinY + cellHalfWidth;

            rc = sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":x"), cellCenterX);
            rc = sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":y"), cellCenterY);

            rc = sqlite3_step(stmt);
            rc = sqlite3_reset(stmt);
        }
    }

    rc = sqlite3_reset(stmt);
    rc = sqlite3_finalize(stmt);

    rc = sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, "PRAGMA SYNCHRONOUS=ON", NULL, NULL, NULL);

    // "Vacuum" the database to free unused memory
    printf("Optimizing database file size, please wait...\n");
    rc = sqlite3_exec(db, "VACUUM", NULL, NULL, NULL);

    rc = sqlite3_close(db);
}
