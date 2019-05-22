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

#include "sqlite3.h"
#include "WfipsGridData.h"

using std::multimap;
using std::string;
using std::vector;
using std::printf;

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

struct MultiPolygonData
{
    OGRSpatialReference spatialReference;
    vector<vector<OGRMultiPolygon>> polygons;
    vector<vector<SBoundingBox>> boundingBoxes;
};

struct WfipsData
{
    OGRSpatialReference spatialReference;
    std::vector<OGRPoint> cellCentroids;
    std::vector<SBoundingBox> cellBoundingBoxes;
    CWfipsGridData gridData;

    double GeoTransform[6];

    int numRows = -1;
    int numCols = -1;
};

struct FireshedData
{
    //vector<int> fireNumbers;
    vector<vector<int>> fireOriginCells;
    vector<unordered_map<int, int>> wfipscellsToFireOriginsForSingleFile;
    vector<vector<int>> originCellsForWfipscellVectorRowIndex;
    unordered_map<int, int> vectorRowIndexToWfipsCell;
    vector<vector<int>> frequenciesForWfipscellVectorRowIndex;
};

bool BoundingBoxCheck(SBoundingBox innerBoundingBox, struct SBoundingBox outerBoundingBox);

void ReadFromShapeFileToMemory(const int shapeFileIndex, const std::string shapeFilePath, const std::string shapeFileName, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData);
void FillFireOriginData(const int shapeFileIndex, string& shapeFileName, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData);
void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData);
int FillWfipsData(WfipsData& wfipsData, std::string dataPath);
void CreateFireShedDB(FireshedData& fireshedData, WfipsData& wfipsData);

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
    int rc = 0;

    if (argc > 1)
    {
        dataPath = argv[1];
    }
    else
    {
        printf("Error: need path to shapefiles as an argument");
         return EXIT_FAILURE;
    }

    if ((dataPath.back() != '/') && (dataPath.back() !='\\'))
    {
        dataPath.push_back('/');
    }

    vector<string> shapeFileNameList;
    vector<string> shapeFilePathList;

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
                    shapeFileNameList.push_back(fileName);
                    shapeFilePathList.push_back(dataPath + fileName);
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

    if (shapeFilePathList.size() < 1)
    {
        printf("Error: No shapefiles found, exiting program\n");
        return EXIT_FAILURE;
    }

    const int num_shape_files = shapeFilePathList.size();

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

    MultiPolygonData multiPolygonData;
    FireshedData fireshedData;
   
    rc = multiPolygonData.spatialReference.importFromESRI(papszPrj5070);
    multiPolygonData.boundingBoxes.resize(num_shape_files);
    multiPolygonData.polygons.resize(num_shape_files);

    fireshedData.fireOriginCells.resize(num_shape_files);
    fireshedData.wfipscellsToFireOriginsForSingleFile.resize(num_shape_files);

    int numCores = omp_get_num_procs();
    if (numCores < 24)
    {
        omp_set_num_threads(numCores - 1);
    }
    else
    {
        numCores = 36;
        omp_set_num_threads(numCores);
    }

    // Parallel
//#pragma omp parallel for
    for (int shapeFileIndex = 0; shapeFileIndex < num_shape_files; shapeFileIndex++)
    {
        ReadFromShapeFileToMemory(shapeFileIndex, dataPath, shapeFilePathList[shapeFileIndex], fireshedData, wfipsData, multiPolygonData);
        FillFireOriginData(shapeFileIndex, shapeFilePathList[shapeFileIndex], fireshedData, wfipsData, multiPolygonData);
    }
    // End Parallel
   
    ConsolidateFinalData(num_shape_files, fireshedData);

    return 0;
}

/***************************************************************************
// Shapefile Reading Function
***************************************************************************/
void ReadFromShapeFileToMemory(const int shapeFileIndex, const std::string shapeFilePath, const std::string shapeFileName, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData)
{
    OGRErr error;
    
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(shapeFileName.c_str(), GDAL_OF_UPDATE, NULL, NULL, NULL));


    OGRLayer  *poLayer = poDS->GetLayer(0);
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();

    OGRwkbGeometryType LayerGeometryType = poLayer->GetGeomType();
    int NumberOfFeatures = poLayer->GetFeatureCount(true);
    poLayer->ResetReading();

    int sizeInAcres = 0;
    int fireNumber;
    SBoundingBox sBoundingBox;
    //Polygon Shapefile

    enum
    {
        fire_number = 0,
        acres = 7,
        x_val = 8,
        y_val = 9
    };

    if (wkbFlatten(LayerGeometryType) == wkbPolygon)
    {
        OGRFeature *poInputFeature;
        
        for (int polygonIndex = 0; polygonIndex < NumberOfFeatures; polygonIndex++)
        {
            poInputFeature = poLayer->GetNextFeature();
            OGRGeometry *poGeometry;

            poGeometry = poInputFeature->GetGeometryRef();
            poGeometry->assignSpatialReference(&(multiPolygonData.spatialReference));
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

            int originCell = 0;
            if (sizeInAcres >= 150)
            {
                originCell = (wfipsData.gridData.WG_GetCellIndex(x, y));
                fireshedData.fireOriginCells[shapeFileIndex].push_back(originCell);
                fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].insert(std::make_pair(originCell, originCell));
                //fireshedData.fireNumbers.push_back(fireNumber);

                if (poGeometry != NULL && (wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon))
                {         
                    OGRMultiPolygon *poPolygon = (OGRMultiPolygon *)poGeometry;
                    poPolygon->closeRings();
                    multiPolygonData.polygons[shapeFileIndex].push_back(*((OGRMultiPolygon*)poPolygon));
                }
                else if (poGeometry != NULL && (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon))
                {
                    poGeometry = OGRGeometryFactory::forceToMultiPolygon(poGeometry);
                    OGRMultiPolygon *poPolygon = (OGRMultiPolygon *)poGeometry;
                    poPolygon->closeRings();
                    multiPolygonData.polygons[shapeFileIndex].push_back(*((OGRMultiPolygon*)poPolygon));
                }

             

                OGREnvelope Envelope;
                poGeometry->getEnvelope(&Envelope);
                sBoundingBox.maxX = Envelope.MaxX;
                sBoundingBox.maxY = Envelope.MaxY;
                sBoundingBox.minX = Envelope.MinX;
                sBoundingBox.minY = Envelope.MinY;
                multiPolygonData.boundingBoxes[shapeFileIndex].push_back(sBoundingBox);
            }
            OGRFeature::DestroyFeature(poInputFeature);
        }
    }

    GDALClose(poDS);
}

bool BoundingBoxCheck(SBoundingBox innerBoundingBox, struct SBoundingBox outerBoundingBox)
{
    return ((innerBoundingBox.minX >= outerBoundingBox.minX) && (innerBoundingBox.maxX <= outerBoundingBox.maxX) && (innerBoundingBox.minY >= outerBoundingBox.minY) && (innerBoundingBox.maxY <= outerBoundingBox.maxY));
}

int FillWfipsData(WfipsData& wfipsData, std::string dataPath)
{
    int rc = wfipsData.gridData.LoadData((const char*)dataPath.c_str());
    std::printf("WFIPS grid data initialized\n");

    //char *MapESRIProjStrings[] =
    //{
    //    "",
    //    "PROJCS[\"Albers\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Albers\"],PARAMETER[\"standard_parallel_1\",29.5],PARAMETER[\"standard_parallel_2\",45.5],PARAMETER[\"latitude_of_origin\",23],PARAMETER[\"central_meridian\",-96],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"METERS\",1]]",
    //    "PROJCS[\"WGS 84 / Pseudo - Mercator\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Mercator\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]"
    //};
    //CPLStringList papszPrj5070 = NULL;
    //OGRSpatialReference spatialReference;
    //papszPrj5070.AddString(MapESRIProjStrings[1]);
    //OGRErr error = spatialReference.importFromESRI(papszPrj5070);

    wfipsData.numRows = wfipsData.gridData.GetNumY();
    wfipsData.numCols = wfipsData.gridData.GetNumX();

    SBoundingBox currentCellBoundingBox;
    int cellIndex = 0;

    wfipsData.gridData.GetWfipsGridGeotransform(wfipsData.GeoTransform);

    for (int row = 0; row < wfipsData.numRows; row++)
    {
      
        for (int col = 0; col < wfipsData.numCols; col++)
        {
            cellIndex = wfipsData.gridData.GetWfipsCellIndex(row, col);
            wfipsData.gridData.WG_GetCellCoords(cellIndex,
                &currentCellBoundingBox.minX,
                &currentCellBoundingBox.minY,
                &currentCellBoundingBox.maxX,
                &currentCellBoundingBox.maxY);
            wfipsData.cellBoundingBoxes.push_back(currentCellBoundingBox);
            OGRPoint cellCenteroid;
            //cellCenteroid.assignSpatialReference(&spatialReference);
            cellCenteroid.setX(currentCellBoundingBox.minX + cellHalfWidth);
            cellCenteroid.setY(currentCellBoundingBox.minY + cellHalfWidth);
            wfipsData.cellCentroids.push_back(cellCenteroid);
        }
    }
    return rc;
}

void FillFireOriginData(const int shapeFileIndex, string& shapeFileName, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData)
{
    int lastError = 0;
    bool isInBoundingBox = false;

    int cellIndex = 0;
    const int num_polygons = multiPolygonData.polygons[shapeFileIndex].size();

    //vector<vector<int>> wfipscellsForEachPolygon;
    //wfipscellsForEachPolygon.resize(num_polygons);

    OGRBoolean intersects = false;

    double currentProgress = 0;
    for (int multipolygonIndex = 0; multipolygonIndex < num_polygons; multipolygonIndex++)
    {
        double minX = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].minX;
        double maxX = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].maxX;
        double minY = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].minY;
        double maxY = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].maxY;

        int upperLeftCell = wfipsData.gridData.WG_GetCellIndex(minX, maxY);
        int lowerRightCell = wfipsData.gridData.WG_GetCellIndex(maxX, minY);

        int minRow = 0;
        int maxRow = 0;
        int minCol = 0;
        int maxCol = 0;

        wfipsData.gridData.GetRowColFromWfipsCell(upperLeftCell, &minRow, &minCol);
        wfipsData.gridData.GetRowColFromWfipsCell(lowerRightCell, &maxRow, &maxCol);

        int origin = fireshedData.fireOriginCells[shapeFileIndex][multipolygonIndex];

        // Create one cell buffer around test if it is within grid bounds
        if (minRow > 1)
        {
            minRow -= 1;
        }
        if (maxRow < wfipsData.numRows)
        {
            maxRow += 1;
        }
        if (minCol > 1)
        {
            minCol -= 1;
        }
        if (maxCol < wfipsData.numCols)
        {
            maxCol += 1;
        }

        for (int row = minRow; row < maxRow; row++)
        {
            for (int col = minCol; col < maxCol; col++)
            {
                cellIndex = wfipsData.gridData.GetWfipsCellIndex(row, col);
                origin = fireshedData.fireOriginCells[shapeFileIndex][multipolygonIndex];
                if (!(cellIndex == origin))
                {
                    isInBoundingBox = BoundingBoxCheck(wfipsData.cellBoundingBoxes[cellIndex], multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex]);
                    if (isInBoundingBox)
                    {
                        SBoundingBox cellBoundingBox = wfipsData.cellBoundingBoxes[cellIndex];
                        OGRLinearRing ring;
                        ring.addPoint(cellBoundingBox.minX, cellBoundingBox.maxY);
                        ring.addPoint(cellBoundingBox.minX, cellBoundingBox.minY);
                        ring.addPoint(cellBoundingBox.maxX, cellBoundingBox.minY);
                        ring.addPoint(cellBoundingBox.maxX, cellBoundingBox.maxY);
                        ring.addPoint(cellBoundingBox.minX, cellBoundingBox.maxY);
                        OGRPolygon cellPoly;
                        cellPoly.addRing(&ring);

                        intersects = cellPoly.Intersects(&multiPolygonData.polygons[shapeFileIndex][multipolygonIndex]);
                        lastError = (int)CPLGetLastErrorNo();

                        if (intersects)
                        {
                            //printf("cell in poly intersection\n");
                            fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].insert(std::make_pair(cellIndex, origin));
                            //wfipscellsForEachPolygon[multipolygonIndex].push_back(cellIndex);
                        }
                    }
                }
            }
        }

        if ((multipolygonIndex > 0) && (multipolygonIndex % 10000 == 0))
        {
            currentProgress = (multipolygonIndex / (num_polygons * 1.0)) * 100.00;
            printf("file\n    %s\n    is %4.2f percent complete\n", shapeFileName.c_str(), currentProgress);
        }
    }

    fireshedData.fireOriginCells[shapeFileIndex].clear();
    //fireshedData.wfipscellsToFireOriginsForSingleFile.clear();

    multiPolygonData.boundingBoxes[shapeFileIndex].clear();
    multiPolygonData.polygons[shapeFileIndex].clear();

    printf("file\n    %s\n   is 100 percent complete\n", shapeFileName.c_str());
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
        for (auto iterator = fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].begin(); iterator != fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].end(); iterator++)
        {
            wfipscell = iterator->first;
            origin = iterator->second;
            totalWfipscellsToFireOrigins.insert(std::make_pair(wfipscell, origin));
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
            fireshedData.vectorRowIndexToWfipsCell.insert(std::make_pair(wfipsCellIndex, wfipscell));
            vector<int> originVectorRow;
            fireshedData.originCellsForWfipscellVectorRowIndex.push_back(originVectorRow);
            vector<int> frequencyVectorRow;
            fireshedData.frequenciesForWfipscellVectorRowIndex.push_back(frequencyVectorRow);
        }
        if (originPrevious != origin) // Origin changed
        {
            originChanged = true;
            originPrevious = origin;
        }

        if (originChanged || wfipsCellChanged) // Add element to current vector row
        {
            frequency = 1;
            fireshedData.originCellsForWfipscellVectorRowIndex[wfipsCellIndex].push_back(origin);
            fireshedData.frequenciesForWfipscellVectorRowIndex[wfipsCellIndex].push_back(frequency);
        }
        else // Increment frequency and overwrite element containing frequency count for current wfipscell
        {
            frequency++;
            fireshedData.frequenciesForWfipscellVectorRowIndex[wfipsCellIndex].back() = frequency;
        }
    }

    totalWfipscellsToFireOrigins.clear();

}

void CreateFireShedDB(FireshedData& fireshedData, WfipsData& wfipsData)
{
    int rc = 0;
    sqlite3 *db;
    sqlite3_stmt *stmt;

    char *sqlErrMsg = 0;

    string createFireshedDBSQLString = "CREATE TABLE IF NOT EXISTS firesheds(" \
        "wfipscell INTEGER," \
        "origin INTEGER,"
        "frequency INTEGER, " \
        "x REAL," \
        "y REAL)";

    string insertTestSQLString = "INSERT INTO firesheds(" \
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
    int origin = -1;
    int err = -1;
    double cellCenterX = -1;
    double cellCenterY = -1;
   
    //for (int wfipsCellIndex = 0; wfipsCellIndex < fireshedData.originCellsForWfipscellVectorRowIndex.size(); wfipsCellIndex++)
    //{
    //    for (int originCellIndex = 0; originCellIndex < fireshedData.originCellsForWfipscellVectorRowIndex[wfipsCellIndex].size(); originCellIndex++)
    //    {
    //        wfipscell = fireshedData.vectorRowIndexToWfipsCell.at(wfipsCellIndex);
    //        bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":wfipscell");
    //        rc = sqlite3_bind_int(stmt, bindColumnIndex, wfipscell);

    //        origin = readOriginCells[wfipsCellIndex][originCellIndex];
    //        bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":origin");
    //        rc = sqlite3_bind_int(stmt, bindColumnIndex, origin);

    //        frequency = readOriginCellFreq[wfipsCellIndex][originCellIndex];
    //        bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":frequency");
    //        rc = sqlite3_bind_int(stmt, bindColumnIndex, frequency);

    //        wfipsData.gridData.WG_GetCellCoords(origin, &cellMinX, &cellMinY, &cellMaxX, &cellMaxY);

    //        double cellCenterX = wfipsData.cellCentroids[wfipsCellIndex].getX();
    //        double cellCenterY = wfipsData.cellCentroids[wfipsCellIndex].getY();

    //        cellCenterX = cellMinX + cellHalfWidth;
    //        cellCenterY = cellMinY + cellHalfWidth;

    //        rc = sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":x"), coords.x);
    //        rc = sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":y"), coords.y);

    //        wgData.WG_GetCellCoords(wfipscell, &cellMinX, &cellMinY, &cellMaxX, &cellMaxY);

    //        cellCenterX = cellMinX + cellHalfWidth;
    //        cellCenterY = cellMinY + cellHalfWidth;

    //        p.setX(cellCenterX);
    //        p.setY(cellCenterY);
    //        p.assignSpatialReference(&source);

    //        err = p.transformTo(&target);

    //        coords.x = p.getX();
    //        coords.y = p.getY();

    //        wfipscellSelfCoords[wfipsCellIndex] = coords;
    //        
    //        rc = sqlite3_step(stmt);
    //        rc = sqlite3_reset(stmt);
    //    }
    //}

    //rc = sqlite3_reset(stmt);
    //rc = sqlite3_finalize(stmt);

    //rc = sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &sqlErrMsg);
    //rc = sqlite3_exec(db, "PRAGMA SYNCHRONOUS=ON", NULL, NULL, NULL);

    //printf("Creating spatial metadata, please wait as this might take a while...\n");
    //rc = sqlite3_exec(db, "SELECT CreateSpatialIndex('test', 'geometry')",
    //    NULL, NULL, &sqlErrMsg);
    //if (rc == SQLITE_OK)
    //{
    //    printf("Spatial metadata created on wfips_cells successfully\n\n");
    //}
    //else
    //{
    //    printf("%s", sqlErrMsg);
    //}

    //// "Vacuum" the database to free unused memory
    //printf("Optimizing database file size, please wait...\n");
    //rc = sqlite3_exec(db, "VACUUM", NULL, NULL, NULL);

    //rc = sqlite3_close(db);
}
