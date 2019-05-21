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

#include "sqlite3.h"
#include "WfipsGridData.h"

using std::multimap;
using std::string;
using std::vector;

/***************************************************************************
// Data Structures
***************************************************************************/

//data structure to hold bounding box
struct SBoundingBox
{
    double MaxX;
    double MaxY;
    double MinX;
    double MinY;
};

struct MultiPolygonData
{
    vector<vector<OGRMultiPolygon>> polygons;
    vector<vector<SBoundingBox>> boundingBoxes;
};

struct WfipsData
{
    std::vector<OGRPoint> cellCentroids;
    std::vector<SBoundingBox> cellBoundingBoxes;
    CWfipsGridData gridData;

    int numRows = -1;
    int numCols = -1;
};

struct FireshedData
{
    //vector<int> fireNumbers;
    vector<vector<int>> fireOriginCells;
    vector<unordered_map<int, int>> wfipscellsToFireOriginsForSingleFile;
    multiset<pair<int, int>> totalWfipscellsToFireOrigins;
    vector<vector<int>> originCellsForWfipscellRowIndex;
    unordered_map<int, int> rowIndexToWfipsCell;
    vector<vector<int>> frequenciesForWfipscellRowIndex;
};

bool BoundingBoxCheck(SBoundingBox innerBoundingBox, struct SBoundingBox outerBoundingBox);

void ReadFromShapeFileToMemory(const int shapeFileIndex, const std::string shapeFilePath, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData);
void FillFireOriginData(const int shapeFileIndex, string& shapeFileName, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData);
void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData);
void FillWfipsData(WfipsData& wfipsData, std::string dataPath);
void CreateFireShedDB(FireshedData& fireshedData, WfipsData& wfipsData);

static const double cellHalfWidth = 1000; // 1 km

int main(int argc, char *argv[])
{
    DIR *dir;
    struct dirent *ent;
    vector<string> shapeFileNameList;
    vector<string> shapeFilePathList;
    string fileName;
    string extension;
    size_t pos;
    int type;

    string dataPath = argv[1];

    if ((dataPath.back() != '/') || (dataPath.back() != '\\'))
    {
        dataPath.push_back('/');
    }

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

    const int num_shape_files = shapeFilePathList.size();

    GDALAllRegister();
   
    WfipsData wfipsData;
    FillWfipsData(wfipsData, dataPath);
    printf("WFIPS data populated\n");

    MultiPolygonData multiPolygonData;
    FireshedData fireshedData;
   
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
        omp_set_num_threads(23);
    }

    // Parallel
#pragma omp parallel for
    for (int shapeFileIndex = 0; shapeFileIndex < num_shape_files; shapeFileIndex++)
    {
        ReadFromShapeFileToMemory(shapeFileIndex, shapeFilePathList[shapeFileIndex], fireshedData, wfipsData, multiPolygonData);
        //printf("polygon and multiPolygon vectors populated for file\n    %s\n", shapeFilePathList[shapeFileIndex].c_str());
        FillFireOriginData(shapeFileIndex, shapeFilePathList[shapeFileIndex], fireshedData, wfipsData, multiPolygonData);
    }
    // End Parallel
   
    ConsolidateFinalData(num_shape_files, fireshedData);

    return 0;
}

/***************************************************************************
// Shapefile Reading Function
***************************************************************************/

void ReadFromShapeFileToMemory(const int shapeFileIndex, const std::string shapeFilePath, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData)
{
    OGRErr error;
    OGRDataSource *poDataSource;

    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(shapeFilePath.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));


    OGRLayer  *poLayer = poDS->GetLayer(0);
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();

    OGRwkbGeometryType LayerGeometryType = poLayer->GetGeomType();
    int NumberOfFeatures = poLayer->GetFeatureCount(true);
    poLayer->ResetReading();

    //char *MapESRIProjStrings[] =
    //{
    //    "",
    //    "PROJCS[\"Albers\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Albers\"],PARAMETER[\"standard_parallel_1\",29.5],PARAMETER[\"standard_parallel_2\",45.5],PARAMETER[\"latitude_of_origin\",23],PARAMETER[\"central_meridian\",-96],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"METERS\",1]]",
    //    "PROJCS[\"WGS 84 / Pseudo - Mercator\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Mercator\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]"
    //};
    //CPLStringList papszPrj5070 = NULL;
    //OGRSpatialReference spatialReference;
    //papszPrj5070.AddString(MapESRIProjStrings[1]);
    //error = spatialReference.importFromESRI(papszPrj5070);

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
        OGRFeature *poFeature;
        for (int polygonIndex = 0; polygonIndex < NumberOfFeatures; polygonIndex++)
        {
            poFeature = poLayer->GetNextFeature();
            OGRGeometry *poGeometry;

            poGeometry = poFeature->GetGeometryRef();

            double x = 0.0;
            double y = 0.0;

            for (int iField = 0; iField < 10; iField++)
            {
                OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);

                if (iField == fire_number)
                {
                    fireNumber = poFeature->GetFieldAsInteger(iField);
                }
                if (iField == acres)
                {
                    sizeInAcres = poFeature->GetFieldAsInteger(iField);
                }
                else if (iField == x_val)
                {
                    x = poFeature->GetFieldAsDouble(iField);
                }
                else if (iField == y_val)
                {
                    y = poFeature->GetFieldAsDouble(iField);
                }
            }

            int originCell = 0;
            if (sizeInAcres >= 150)
            {
                originCell = (wfipsData.gridData.WG_GetCellIndex(x, y));
                fireshedData.fireOriginCells[shapeFileIndex].push_back(originCell);
                fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].insert(std::make_pair(originCell, originCell));
                //multiPolygonData.fireNumbers.push_back(fireNumber);

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
                sBoundingBox.MaxX = Envelope.MaxX;
                sBoundingBox.MaxY = Envelope.MaxY;
                sBoundingBox.MinX = Envelope.MinX;
                sBoundingBox.MinY = Envelope.MinY;
                multiPolygonData.boundingBoxes[shapeFileIndex].push_back(sBoundingBox);
            }
            OGRFeature::DestroyFeature(poFeature);
        }
    }

    GDALClose(poDS);
}

bool BoundingBoxCheck(SBoundingBox innerBoundingBox, struct SBoundingBox outerBoundingBox)
{
    return ((innerBoundingBox.MinX >= outerBoundingBox.MinX) && (innerBoundingBox.MaxX <= outerBoundingBox.MaxX) && (innerBoundingBox.MinY >= outerBoundingBox.MinY) && (innerBoundingBox.MaxY <= outerBoundingBox.MaxY));
}

void FillWfipsData(WfipsData& wfipsData, std::string dataPath)
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

    int cellIndex = 0;
    for (int row = 0; row < wfipsData.numRows; row++)
    {
        SBoundingBox currentCellBoundingBox;
        for (int col = 0; col < wfipsData.numCols; col++)
        {
            cellIndex = wfipsData.gridData.GetWfipsCellIndex(row, col);
            wfipsData.gridData.WG_GetCellCoords(cellIndex,
                &currentCellBoundingBox.MinX,
                &currentCellBoundingBox.MinY,
                &currentCellBoundingBox.MaxX,
                &currentCellBoundingBox.MaxY);
            wfipsData.cellBoundingBoxes.push_back(currentCellBoundingBox);
            OGRPoint cellCenteroid;
            //cellCenteroid.assignSpatialReference(&spatialReference);
            cellCenteroid.setX(currentCellBoundingBox.MinX + cellHalfWidth);
            cellCenteroid.setY(currentCellBoundingBox.MinY + cellHalfWidth);
            wfipsData.cellCentroids.push_back(cellCenteroid);
        }
    }
}

void FillFireOriginData(const int shapeFileIndex, string& shapeFileName, FireshedData& fireshedData, WfipsData& wfipsData, MultiPolygonData& multiPolygonData)
{
    int lastError = 0;
    bool isInBoundingBox = false;
    OGRBoolean isIntersecting = false;

    int cellIndex = 0;
    const int num_polygons = multiPolygonData.polygons[shapeFileIndex].size();

    double currentProgress = 0;
    for (int multipolygonIndex = 0; multipolygonIndex < num_polygons; multipolygonIndex++)
    {
        double minX = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].MinX;
        double maxX = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].MaxX;
        double minY = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].MinY;
        double maxY = multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex].MaxY;

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

                if (!(cellIndex == fireshedData.fireOriginCells[shapeFileIndex][multipolygonIndex]))
                {
                    isInBoundingBox = BoundingBoxCheck(wfipsData.cellBoundingBoxes[cellIndex], multiPolygonData.boundingBoxes[shapeFileIndex][multipolygonIndex]);
                    if (isInBoundingBox)
                    {
                        isIntersecting = wfipsData.cellCentroids[cellIndex].Intersects(&multiPolygonData.polygons[shapeFileIndex][multipolygonIndex]);
                        lastError = (int)CPLGetLastErrorNo();

                        if (isIntersecting)
                        {
                            //printf("Intersection found\n");
                            fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].insert(std::make_pair(cellIndex, origin));
                        }
                    }
                }
            }
        }

        if ((multipolygonIndex > 0) && (multipolygonIndex % 10000 == 0))
        {
            currentProgress = (multipolygonIndex / (num_polygons * 1.0)) * 100.00;
            printf("file %s\n    is %4.2f percent complete\n", shapeFileName.c_str(), currentProgress);
        }
    }

    fireshedData.fireOriginCells[shapeFileIndex].clear();
    //fireshedData.wfipscellsToFireOriginsForSingleFile.clear();

    multiPolygonData.boundingBoxes[shapeFileIndex].clear();
    multiPolygonData.polygons[shapeFileIndex].clear();

    printf("file %s\n    is 100 percent complete\n", shapeFileName.c_str(), currentProgress);
}

void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData)
{
    multiset<pair<int, int>> totalWfipscellsToFireOrigins;

    vector<vector<int>> originCellsForWfipscellRowIndex;
    unordered_map<int, int> rowIndexToWfipsCell;
    vector<vector<int>> frequenciesForWfipscellRowIndex;

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
            rowIndexToWfipsCell.insert(std::make_pair(wfipsCellIndex, wfipscell));
            vector<int> originVectorRow;
            originCellsForWfipscellRowIndex.push_back(originVectorRow);
            vector<int> frequencyVectorRow;
            frequenciesForWfipscellRowIndex.push_back(frequencyVectorRow);
        }
        if (originPrevious != origin) // Origin changed
        {
            originChanged = true;
            originPrevious = origin;
        }

        if (originChanged || wfipsCellChanged) // Add element to current vector row
        {
            frequency = 1;
            originCellsForWfipscellRowIndex[wfipsCellIndex].push_back(origin);
            frequenciesForWfipscellRowIndex[wfipsCellIndex].push_back(frequency);
        }
        else // Increment frequency and overwrite element containing frequency count for current wfipscell
        {
            frequency++;
            frequenciesForWfipscellRowIndex[wfipsCellIndex].back() = frequency;
        }
    }

    totalWfipscellsToFireOrigins.clear();

}

void CreateFireShedDB(FireshedData& fireshedData, WfipsData& wfipsData)
{
    //for (int wfipsCellIndex = 0; wfipsCellIndex < readOriginCells.size(); wfipsCellIndex++)
    //{
    //    for (int originCellIndex = 0; originCellIndex < readOriginCells[wfipsCellIndex].size(); originCellIndex++)
    //    {
    //        wfipscell = rowIndexToWfipsCell.at(wfipsCellIndex);
    //        bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":wfipscell");
    //        rc = sqlite3_bind_int(stmt, bindColumnIndex, wfipscell);

    //        origin = readOriginCells[wfipsCellIndex][originCellIndex];
    //        bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":origin");
    //        rc = sqlite3_bind_int(stmt, bindColumnIndex, origin);

    //        frequency = readOriginCellFreq[wfipsCellIndex][originCellIndex];
    //        bindColumnIndex = sqlite3_bind_parameter_index(stmt, ":frequency");
    //        rc = sqlite3_bind_int(stmt, bindColumnIndex, frequency);

    //        wfipsData.gridData.WG_GetCellCoords(origin, &cellMinX, &cellMinY, &cellMaxX, &cellMaxY);

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
