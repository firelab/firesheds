
#include <iostream>

#ifdef  _WIN32
#include "win\dirent.h"
#else
#include "dirent.h"
#endif

#include "WfipsGrid.h"
#include "ogrsf_frmts.h"
#include "gdal_alg.h"

#include <algorithm>
#include <omp.h>
#include <set>
#include <vector>
#include <unordered_map>
#include <cstdio>
#include <chrono>
#include <gdal_utils.h>

#include "sqlite3.h"

#include "WfipsGridData.h"

#include "MyPolygon.h"

#ifndef EQUAL
#define EQUAL(a,b) (strcmp(a,b)==0)
#endif

using std::multimap;
using std::string;
using std::vector;
using std::chrono::steady_clock;

/***************************************************************************
// Global Constants
***************************************************************************/
static const double epsilon = 0.000000000000000001;
static const double cellWidthInMeters = 2000.0; // 2 km
static const double numMicrosecondsPerSecond = 1000000.0;

/***************************************************************************
// Data Structures
***************************************************************************/

struct WfipsData
{
    OGRSpatialReference spatialReference;
    vector<SBoundingBox> cellBoundingBoxes;
    CWfipsGridData gridData;

    double GeoTransform[6];

    int numRows = -1;
    int numCols = -1;
};

struct FireshedData
{
    vector<unordered_multimap<int, int>> wfipscellsToFireOriginsForSingleFile;
    vector<vector<int>> originCellsForWfipscell;
    unordered_map<int, int> finalIndexToWfipsCellMap;
    vector<vector<int>> numWfipscellOriginPairs;
    map<int, int> finalOriginCellToTotalPairCountMap;
};

SBoundingBox GetBoundingBox(vector<MyPoint2D> theRingString);

void ReadShapefilesToMemory(const int numThreads, const bool verbose, const steady_clock::time_point startClock, const int numEdgeCellDivisions, string& shapefilePath, vector<string>& shapefileNameList, FireshedData& fireshedData, WfipsData& wfipsData);
void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData);
int FillWfipsData(WfipsData& wfipsData, std::string dataPath);
int CreateFireShedDB(const bool verbose, sqlite3* db, FireshedData& fireshedData, WfipsData& wfipsData);
void Usage();

bool AreClose(double a, double b);



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
    string gridPath = "";
    string numThreadsString = "";
    int rc = 0;
    bool verboseParameter = false;
    const int SUCCESS = 0;

    const int max_argument_index = argc - 1;
    int argIndex = 1;
    bool isOutPathSpecified = false;
    bool isGridPathSpecified = false;
    bool isNumThreadsSpecified = false;

    // Parse commandline arguments
    if (argc > 1)
    {
        while (argIndex < argc)
        {
            if (EQUAL(argv[argIndex], "--shape"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No shapefile name entered\n");
                    Usage(); // Exits program
                }
                dataPath = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--out"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No output path entered\n");
                    Usage(); // Exits program
                }
                isOutPathSpecified = true;
                outPath = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--grid"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No grid path entered\n");
                    Usage(); // Exits program
                }
                isGridPathSpecified = true;
                gridPath = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--t"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No thread number entered\n");
                    Usage(); // Exits program
                }
                isNumThreadsSpecified = true;
                numThreadsString = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--verbose"))
            {
                verboseParameter = true;
            }
            else
            {
                printf("ERROR: %s is an invalid argument\n", argv[argIndex]);
                Usage(); // Exits program
            }
            argIndex++;
        }
    }
    else
    {
        printf("ERROR: no shapefile path given\n", argv[argIndex]);
        Usage(); // Exits program
    }

    const bool verbose = verboseParameter;

    if ((dataPath.back() != '/') && (dataPath.back() != '\\'))
    {
#ifdef WIN32
        dataPath.push_back('\\');
#else
        dataPath.push_back('/');
#endif
    }

    if (!isGridPathSpecified)
    {
        gridPath = dataPath;
    }
    else if ((gridPath.back() != '/') && (gridPath.back() != '\\'))
    {
#ifdef WIN32
        gridPath.push_back('\\');
#else
        gridPath.push_back('/');
#endif
    }

    if (!isOutPathSpecified)
    {
        outPath = dataPath;
    }
    else if ((outPath.back() != '/') && (outPath.back() != '\\'))
    {
#ifdef WIN32
        outPath.push_back('\\');
#else
        outPath.push_back('/');
#endif
    }

    int numThreadsArg = -1;
    if (isNumThreadsSpecified)
    {
        numThreadsArg = atoi(numThreadsString.c_str());
    }

    if (isNumThreadsSpecified && (numThreadsArg < 1))
    {
        printf("ERROR: Number of threads must be greater than zero\n", argv[argIndex]);
        Usage(); // Exits program
    }

    vector<string> shapefileNameList;
    vector<string> shapefilePathList;

    if ((dir = opendir(dataPath.c_str())) != nullptr)
    {
        /* readt all the files and directories within directory */
        while ((ent = readdir(dir)) != nullptr)
        {
            fileName = ent->d_name;

            pos = fileName.find_last_of('.');
            if (pos != string::npos)
            {
                extension = fileName.substr(pos, fileName.size());
                std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
                if (extension == ".shp")
                {
                    //printf("%s\n", fileName.c_str());
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

    sqlite3 *db = nullptr;

    string outPutDbFullPath = outPath + "firesheds.db";
    rc = sqlite3_open_v2(outPutDbFullPath.c_str(), &db,
        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
        NULL);

    if (rc != SQLITE_OK || db == nullptr)
    {
        printf("Error: Could not create firesheds.db (does the output directory exist?)\nexiting program\n");
        return EXIT_FAILURE;
    }

    GDALAllRegister();

    char *MapESRIProjStrings[] =
    {
        "",
        "PROJCS[\"Albers\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Albers\"],PARAMETER[\"standard_parallel_1\",29.5],PARAMETER[\"standard_parallel_2\",45.5],PARAMETER[\"latitude_of_origin\",23],PARAMETER[\"central_meridian\",-96],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"METERS\",1]]"
    };

    CPLStringList papszPrj5070;
    papszPrj5070.AddString(MapESRIProjStrings[1]);

    WfipsData wfipsData;
    rc = wfipsData.spatialReference.importFromESRI(papszPrj5070);

    if (verbose)
    {
        printf("Loading WFIPS data\n");
    }

    rc = FillWfipsData(wfipsData, dataPath);
    if (rc != SUCCESS)
    {
        printf("Error: WFIPS data loading failed\n");
        printf("    Make sure the file \"WFIPSGrid.tif\" exists in\n    %s\n\n", dataPath.c_str());
        return EXIT_FAILURE;
    }
    else if (verbose)
    {
        printf("WFIPS data populated\n\n");
    }

    printf("Processing all shapefiles in\n    %s\n\nPlease wait...\n\n", dataPath.c_str());

    FireshedData fireshedData;
    fireshedData.wfipscellsToFireOriginsForSingleFile.resize(shapefileListSize);

    steady_clock::time_point startClock = std::chrono::steady_clock::now();
    const int numCellEdgeDivisions = 16; // Number of times to subdivide WFIPS cells (currently 2000m, makes 125m subcells) 
                                         // to get closer to the resolution used in FSim (apparently ~135m)

    ReadShapefilesToMemory(numThreadsArg, verbose, startClock, numCellEdgeDivisions, dataPath, shapefileNameList, fireshedData, wfipsData);

    ConsolidateFinalData(shapefileListSize, fireshedData);
    CreateFireShedDB(verbose, db, fireshedData, wfipsData);

    long long int elapsedTimeInMicroseconds = std::chrono::microseconds(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startClock)).count();
    printf("Successfully processed all shapefiles in\n    %s\n", dataPath.c_str());
    printf("    and created \"firesheds.db\" in\n    %s\n\n", outPath.c_str());
    printf("Total time elapsed is %4.2f seconds\n\n", elapsedTimeInMicroseconds / numMicrosecondsPerSecond);

    return SUCCESS;
}

/***************************************************************************
// Shapefile Reading Function
***************************************************************************/
void ReadShapefilesToMemory(const int numThreadsArg, const bool verbose, const steady_clock::time_point startClock, const int numEdgeCellDivisions, string& shapefilePath, vector<string>& shapefileNameList, FireshedData& fireshedData, WfipsData& wfipsData)
{
    int numThreads = 1;
    const int availableThreads = omp_get_num_procs();
    if (numThreadsArg == -1)
    {
        numThreads = availableThreads;
    }
    else
    {
        if (numThreadsArg > availableThreads)
        {
            numThreads = availableThreads;
        }
        else
        {
            numThreads = numThreadsArg;
        }
    }

    omp_set_num_threads(numThreads);

    int numFilesProcessed = 0;

    const int shapefileListSize = shapefileNameList.size();

    #pragma omp parallel for schedule(dynamic, 1) shared(numFilesProcessed)
    for (int shapefileIndex = 0; shapefileIndex < shapefileListSize; shapefileIndex++)
    {
        string shapefileFullPath = shapefilePath + shapefileNameList[shapefileIndex];
        GDALDataset *poShapefileDS = static_cast<GDALDataset*>(GDALOpenEx(shapefileFullPath.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL));

        OGRLayer  *poLayer = poShapefileDS->GetLayer(0);
        OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();

        OGRwkbGeometryType LayerGeometryType = poLayer->GetGeomType();
        int NumberOfFeatures = poLayer->GetFeatureCount(true);

        //vector<int> fireNumbers;
        vector<int> fireOriginCells;

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

        //Bounding Box of Shapefile 
        SBoundingBox sBoundingBox;

        //Holds Coordinates of Polygon Shapefile
        vector<PolygonFeature> PolygonLayer;

        int numFields = 10;
        double x = 0.0;
        double y = 0.0;

        int originCell = 0;

        //Polygon Shapefile
        if (wkbFlatten(LayerGeometryType) == wkbPolygon)
        {
            OGRFeature *poFeature;
            PolygonFeature Polygon;
            OGRPoint ptTemp;
            for (int i = 0; i < NumberOfFeatures; i++)
            {
                poFeature = poLayer->GetNextFeature();
                OGRGeometry *poGeometry;
                poGeometry = poFeature->GetGeometryRef();
                poGeometry->assignSpatialReference(&wfipsData.spatialReference);

                for (int fieldIndex = 0; fieldIndex < numFields; fieldIndex++)
                {
                    OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(fieldIndex);

                    if (fieldIndex == fire_number)
                    {
                        fireNumber = poFeature->GetFieldAsInteger(fieldIndex);
                    }
                    if (fieldIndex == acres)
                    {
                        sizeInAcres = poFeature->GetFieldAsInteger(fieldIndex);
                    }
                    else if (fieldIndex == x_val)
                    {
                        x = poFeature->GetFieldAsDouble(fieldIndex);
                    }
                    else if (fieldIndex == y_val)
                    {
                        y = poFeature->GetFieldAsDouble(fieldIndex);
                    }
                }

                if (sizeInAcres >= 150)
                {
                    originCell = wfipsData.gridData.WG_GetCellIndex(x, y);
                    //fireNumbers.push_back(fireNumber);
                    fireOriginCells.push_back(originCell);

                    if (poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon)
                    {
                        OGRPolygon *poPolygon = (OGRPolygon *)poGeometry;
                        Polygon.PolygonsOfFeature.resize(1);
                        int NumberOfInnerRings = poPolygon->getNumInteriorRings();
                        OGRLinearRing *poExteriorRing = poPolygon->getExteriorRing();
                        Polygon.PolygonsOfFeature.at(0).Polygon.resize(NumberOfInnerRings + 1);
                        Polygon.PolygonsOfFeature.at(0).Polygon.at(0).IsClockwised = poExteriorRing->isClockwise();
                        int NumberOfExteriorRingVertices = poExteriorRing->getNumPoints();
                        Polygon.PolygonsOfFeature.at(0).Polygon.at(0).RingString.resize(NumberOfExteriorRingVertices);
                        for (int k = 0; k < NumberOfExteriorRingVertices; k++)
                        {
                            poExteriorRing->getPoint(k, &ptTemp);
                            MyPoint2D pt;
                            pt.X = ptTemp.getX();
                            pt.Y = ptTemp.getY();
                            Polygon.PolygonsOfFeature.at(0).Polygon.at(0).RingString.at(k) = pt;
                        }
                        for (int h = 1; h <= NumberOfInnerRings; h++)
                        {
                            OGRLinearRing *poInteriorRing = poPolygon->getInteriorRing(h - 1);
                            Polygon.PolygonsOfFeature.at(0).Polygon.at(h).IsClockwised = poInteriorRing->isClockwise();
                            int NumberOfInteriorRingVertices = poInteriorRing->getNumPoints();
                            Polygon.PolygonsOfFeature.at(0).Polygon.at(h).RingString.resize(NumberOfInteriorRingVertices);
                            for (int k = 0; k < NumberOfInteriorRingVertices; k++)
                            {
                                poInteriorRing->getPoint(k, &ptTemp);
                                MyPoint2D pt;
                                pt.X = ptTemp.getX();
                                pt.Y = ptTemp.getY();
                                Polygon.PolygonsOfFeature.at(0).Polygon.at(h).RingString.at(k) = pt;
                            }
                        }
                        PolygonLayer.push_back(Polygon);
                    }
                    else if (poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon)
                    {
                        OGRMultiPolygon *poMultiPolygon = (OGRMultiPolygon *)poGeometry;
                        int NumberOfGeometries = poMultiPolygon->getNumGeometries();
                        Polygon.PolygonsOfFeature.resize(NumberOfGeometries);
                        for (int j = 0; j < NumberOfGeometries; j++)
                        {
                            OGRGeometry *poPolygonGeometry = poMultiPolygon->getGeometryRef(j);
                            OGRPolygon *poPolygon = (OGRPolygon *)poPolygonGeometry;
                            int NumberOfInnerRings = poPolygon->getNumInteriorRings();
                            OGRLinearRing *poExteriorRing = poPolygon->getExteriorRing();
                            Polygon.PolygonsOfFeature.at(j).Polygon.resize(NumberOfInnerRings + 1);
                            Polygon.PolygonsOfFeature.at(j).Polygon.at(0).IsClockwised = poExteriorRing->isClockwise();
                            int NumberOfExteriorRingVertices = poExteriorRing->getNumPoints();
                            Polygon.PolygonsOfFeature.at(j).Polygon.at(0).RingString.resize(NumberOfExteriorRingVertices);
                            for (int k = 0; k < NumberOfExteriorRingVertices; k++)
                            {
                                poExteriorRing->getPoint(k, &ptTemp);
                                MyPoint2D pt;
                                pt.X = ptTemp.getX();
                                pt.Y = ptTemp.getY();
                                Polygon.PolygonsOfFeature.at(j).Polygon.at(0).RingString.at(k) = pt;
                            }
                            for (int h = 1; h <= NumberOfInnerRings; h++)
                            {
                                OGRLinearRing *poInteriorRing = poPolygon->getInteriorRing(h - 1);
                                Polygon.PolygonsOfFeature.at(j).Polygon.at(h).IsClockwised = poInteriorRing->isClockwise();
                                int NumberOfInteriorRingVertices = poInteriorRing->getNumPoints();
                                Polygon.PolygonsOfFeature.at(j).Polygon.at(h).RingString.resize(NumberOfInteriorRingVertices);
                                for (int k = 0; k < NumberOfInteriorRingVertices; k++)
                                {
                                    poInteriorRing->getPoint(k, &ptTemp);
                                    MyPoint2D pt;
                                    pt.X = ptTemp.getX();
                                    pt.Y = ptTemp.getY();
                                    Polygon.PolygonsOfFeature.at(j).Polygon.at(h).RingString.at(k) = pt;
                                }
                            }
                        }
                        PolygonLayer.push_back(Polygon);
                    }
                }
                OGRFeature::DestroyFeature(poFeature);
            }
        }

        //fireNumbers.shrink_to_fit();
        fireOriginCells.shrink_to_fit();

        GDALClose(poShapefileDS);

        int numFires = PolygonLayer.size();

        unordered_multimap<int, int> tempTotalWfipscellsToFireOrigins;

        for (int fireIndex = 0; fireIndex < PolygonLayer.size(); fireIndex++)
        {
            //fireNumber = fireNumbers[fireIndex];
            int origin = fireOriginCells[fireIndex];

            unordered_map<int, int> tempSingleFireWfipscellsToFireOrigins; // used to check if cell is already burned by fire

             // Always assume cell containing the origin is burned by the fire
            tempTotalWfipscellsToFireOrigins.insert(std::make_pair(origin, origin));
            tempSingleFireWfipscellsToFireOrigins.insert(std::make_pair(origin, origin));

            int numMultiPolygons = PolygonLayer[fireIndex].PolygonsOfFeature.size();
            for (int multiPolygonIndex = 0; multiPolygonIndex < numMultiPolygons; multiPolygonIndex++)
            {
                int numPolygons = PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon.size();
                for (int polygonIndex = 0; polygonIndex < numPolygons; polygonIndex++)
                {
                    int ringStringSize = PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString.size();
                    MyPoint2D firstPt;
                    firstPt.X = PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString[0].X;
                    firstPt.Y = PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString[0].Y;
                    MyPoint2D lastPt;
                    lastPt.X = PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString[ringStringSize - 1].X;
                    lastPt.Y = PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString[ringStringSize - 1].Y;

                    if (!((AreClose(firstPt.X, lastPt.X)) && (AreClose(firstPt.Y, lastPt.Y))))
                    {
                        PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString.push_back(firstPt);
                    }

                    // Get the bounding box for the current fire
                    SBoundingBox fireBoundingBox = GetBoundingBox(PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString);
                    int upperLeftCell = wfipsData.gridData.WG_GetCellIndex(fireBoundingBox.minX, fireBoundingBox.maxY);
                    int lowerRightCell = wfipsData.gridData.WG_GetCellIndex(fireBoundingBox.maxX, fireBoundingBox.minY);

                    int minRow = 0;
                    int maxRow = 0;
                    int minCol = 0;
                    int maxCol = 0;

                    // Find rows and columns that correspond to the upper-left and lower-right corners of fire bounding box
                    wfipsData.gridData.GetRowColFromWfipsCell(upperLeftCell, &minRow, &minCol);
                    wfipsData.gridData.GetRowColFromWfipsCell(lowerRightCell, &maxRow, &maxCol);

                    int cellIndex = 0;
                    for (int row = minRow; row <= maxRow; row++)
                    {
                        for (int col = minCol; col <= maxCol; col++)
                        {
                            cellIndex = wfipsData.gridData.GetWfipsCellIndex(row, col);
                            origin = fireOriginCells[fireIndex];

                            // Below each cell within the current fire's bounding box not already marked as burned by that fire
                            // is subdivided into numEdgeCellDivisions subcells. The nodes along the outer edges of the subcells
                            // are then used as test points in Mark's point in poly function IsOverlapping(). We test only points
                            // along the edge as we only need to check where the fire crosses into a WFIPS cell from another cell
                            // as origin cells are always assumed to be burned by a fire. If any node on the edge of the subcells
                            // for the current WFIPS cell is found to be within the fire, we mark that cell as burned and move on
                            // to the next cell within the fire's bounding box.

                            // Need to check if cell has already been burned by one of the current fire's subpolygons so as not to
                            // count cells as being burned by the same fire more than once  
                            bool isCellAlreadyBurned = !(tempSingleFireWfipscellsToFireOrigins.find(cellIndex) == tempSingleFireWfipscellsToFireOrigins.end());
                            if (!isCellAlreadyBurned)
                            {
                                SBoundingBox cellBoundingBox = wfipsData.cellBoundingBoxes[cellIndex];

                                const double xMin = cellBoundingBox.minX;
                                const double yMin = cellBoundingBox.minY;
                                const double xMax = cellBoundingBox.maxX;
                                const double yMax = cellBoundingBox.maxY;

                                MyPoint2D testPoint;
                                testPoint.X = xMin;
                                testPoint.Y = yMin;
                                const int numEdgeNodes = numEdgeCellDivisions + 1;
                                const double nodeStepSizeInMeters = cellWidthInMeters * (1.0 / numEdgeCellDivisions);

                                bool isCellInFirePolygon = false;
                                double xEdgeNodeCoord = 0;
                                double yEdgeNodeCoord = 0;

                                const int topEdge = 0; // yNodeIndex value for top edge of cell
                                const int bottomEdge = numEdgeCellDivisions; // yNodeIndex value for bottom edge of WFIPS cell
                                const int leftEdge = 0; // xNodeIndex value for left edge of cell
                                const int rightEdge = numEdgeCellDivisions; // xNodeIndex value for right edge of WFIPS cell

                                for (int xNodeIndex = 0; xNodeIndex < numEdgeNodes; xNodeIndex++)
                                {
                                    if (!isCellInFirePolygon)
                                    {
                                        xEdgeNodeCoord = xMin + (nodeStepSizeInMeters * xNodeIndex);
                                        testPoint.X = xEdgeNodeCoord;
                                        for (int yNodeIndex = 0; yNodeIndex < numEdgeNodes; yNodeIndex++)
                                        {
                                            yEdgeNodeCoord = yMin + (nodeStepSizeInMeters * yNodeIndex);
                                            testPoint.Y = yEdgeNodeCoord;

                                            if ((yNodeIndex == topEdge) || (yNodeIndex == bottomEdge))
                                            {
                                                // Check all subcell nodes along the cell top or bottom egdes of WFIPS cell  
                                                isCellInFirePolygon = MyPolygonUtility::IsOverlapping(testPoint,
                                                    PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString);
                                                if (isCellInFirePolygon)
                                                {
                                                    break; // Leave yNode loop
                                                }
                                            }
                                            else 
                                            {
                                                // If not along top or bottom edge of WFIPS cell, check only subcell nodes on left and right edges
                                                if ((xNodeIndex == leftEdge) || (xNodeIndex == rightEdge))
                                                {
                                                    isCellInFirePolygon = MyPolygonUtility::IsOverlapping(testPoint,
                                                        PolygonLayer[fireIndex].PolygonsOfFeature[multiPolygonIndex].Polygon[polygonIndex].RingString);
                                                    if (isCellInFirePolygon)
                                                    {
                                                        break; // Leave yNode loop
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        break;  // Leave xNode loop
                                    }
                                }
                                if (isCellInFirePolygon)
                                {
                                    // Mark current cell as being burned by the fire
                                    tempTotalWfipscellsToFireOrigins.insert(std::make_pair(cellIndex, origin));
                                    tempSingleFireWfipscellsToFireOrigins.insert(std::make_pair(cellIndex, origin));
                                }
                            }
                        }
                    }
                }
            }
            // Clear list of cells burned by current fire for next loop iteration
            tempSingleFireWfipscellsToFireOrigins.clear();
        }

        int wfipscell = -1;
        int origin = -1;

        #pragma omp atomic
        numFilesProcessed++;

        #pragma omp critical
        {
            // Populate data in shared vector
            for (auto iterator = tempTotalWfipscellsToFireOrigins.begin(); iterator != tempTotalWfipscellsToFireOrigins.end(); iterator++)
            {
                wfipscell = iterator->first;
                origin = iterator->second;
                fireshedData.wfipscellsToFireOriginsForSingleFile[shapefileIndex].insert(std::make_pair(wfipscell, origin));
            }

            if (verbose)
            {
                double currentProgress = (numFilesProcessed / (shapefileListSize * 1.0)) * 100.00;
                long long int elapsedTimeInMicroseconds = std::chrono::microseconds(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startClock)).count();
                printf("Processed %d files out of %d in\n", numFilesProcessed, shapefileListSize);
                printf("    %s\n    %4.2f percent of all files to be processed are complete\n", shapefilePath.c_str(), currentProgress);
                printf(" Total time elapsed is %4.2f seconds\n\n", elapsedTimeInMicroseconds / numMicrosecondsPerSecond);
            }
        }

        //fireNumbers.clear();
        fireOriginCells.clear();
        tempTotalWfipscellsToFireOrigins.clear();
    }
}

int FillWfipsData(WfipsData& wfipsData, std::string dataPath)
{
    int rc = wfipsData.gridData.LoadData((const char*)dataPath.c_str());

    wfipsData.numRows = wfipsData.gridData.GetNumY();
    wfipsData.numCols = wfipsData.gridData.GetNumX();

    wfipsData.gridData.GetWfipsGridGeotransform(wfipsData.GeoTransform);

    SBoundingBox currentCellBoundingBox;

    int cellIndex = 0;
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
        }
    }

    return rc;
}

void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData)
{
    multiset<pair<int, int>> totalWfipscellsToFireOriginPairs;
    vector<pair<int, int>> totalPairsForSingleOrigin;
    int wfipscell = -1;
    int origin = -1;
    int wfipscellPrevious = -1;
    int originPrevious = -1;
    int wfipsCellIndex = -1;
    int numPairs = -1;
    int numTotalPairs = -1;

    for (int shapeFileIndex = 0; shapeFileIndex < num_shape_files; shapeFileIndex++)
    {
        for (auto iterator = fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].begin(); iterator != fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].end(); iterator++)
        {
            wfipscell = iterator->first;
            origin = iterator->second;
            totalWfipscellsToFireOriginPairs.insert(std::make_pair(wfipscell, origin));
        }

        fireshedData.wfipscellsToFireOriginsForSingleFile[shapeFileIndex].clear();
    }

    // Get the data for num_pairs field
    for (auto iterator = totalWfipscellsToFireOriginPairs.begin(); iterator != totalWfipscellsToFireOriginPairs.end(); iterator++)
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
            fireshedData.finalIndexToWfipsCellMap.insert(std::make_pair(wfipsCellIndex, wfipscell));
            vector<int> originVectorRow;
            fireshedData.originCellsForWfipscell.push_back(originVectorRow);
            vector<int> numPairsVectorRow;
            fireshedData.numWfipscellOriginPairs.push_back(numPairsVectorRow);
        }
        if (origin != originPrevious) // Origin changed
        {
            originChanged = true;
            originPrevious = origin;
        }

        if (originChanged || wfipsCellChanged) // Add element to current vector row
        {
            numPairs = 1;
            fireshedData.originCellsForWfipscell[wfipsCellIndex].push_back(origin);
            fireshedData.numWfipscellOriginPairs[wfipsCellIndex].push_back(numPairs);
        }
        else // Increment count of pairs and overwrite element containing count for current wfipscell
        {
            numPairs++;
            fireshedData.numWfipscellOriginPairs[wfipsCellIndex].back() = numPairs;
        }
    }

    // Get the data for total_for_origin field
    for (int wfipsCellIndex = 0; wfipsCellIndex < fireshedData.finalIndexToWfipsCellMap.size(); wfipsCellIndex++)
    {
        int wfipscell = fireshedData.finalIndexToWfipsCellMap.at(wfipsCellIndex);
        for (int originCellIndex = 0; originCellIndex < fireshedData.originCellsForWfipscell[wfipsCellIndex].size(); originCellIndex++)
        {
            origin = fireshedData.originCellsForWfipscell[wfipsCellIndex][originCellIndex];
            if (origin == wfipscell) // Fire started in that cell
            {
                int numTotalPairs = fireshedData.numWfipscellOriginPairs[wfipsCellIndex][originCellIndex];
                fireshedData.finalOriginCellToTotalPairCountMap.insert(std::make_pair(origin, numTotalPairs));
            }
        }
    }

    totalWfipscellsToFireOriginPairs.clear();
}

int CreateFireShedDB(const bool verbose, sqlite3* db, FireshedData& fireshedData, WfipsData& wfipsData)
{
    int rc = 0;
    sqlite3_stmt *fireshedsStmt;
    sqlite3_stmt *totalsStmt;
    char *sqlErrMsg = 0;

    string createFireshedDBSQLString = "CREATE TABLE IF NOT EXISTS firesheds " \
        "(wfipscell INTEGER, " \
        "origin INTEGER, " \
        "num_pairs INTEGER)";
    
    string insertFireshedsSQLString = "INSERT INTO firesheds(" \
        "wfipscell, " \
        "origin, " \
        "num_pairs) " \
        "VALUES " \
        "(:wfipscell, " \
        ":origin, " \
        ":num_pairs)";

    string createTotalFiresDBSQLString = "CREATE TABLE IF NOT EXISTS totals_for_origins " \
        "(origin INTEGER, " \
        "total_fires INTEGER)";

    string insertTotalsFiresSQLString = "INSERT INTO totals_for_origins " \
        "(origin, " \
        "total_fires) " \
        "VALUES " \
        "(:origin, " \
        ":total_fires)";

    int wfipscell = -1;
    int bindColumnIndex = -1;
    int numPairs = -1;
    int origin = -1;


    rc = sqlite3_exec(db, "DROP TABLE IF EXISTS 'firesheds'", NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, "DROP TABLE IF EXISTS 'totals_for_origins'", NULL, NULL, &sqlErrMsg);

    rc = sqlite3_exec(db, createFireshedDBSQLString.c_str(), NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, createTotalFiresDBSQLString.c_str(), NULL, NULL, &sqlErrMsg);

    if (rc != SQLITE_OK)
    {
        printf("Error: Could not create schema on firesheds.db\n\n");
        return rc;
    }

    rc = sqlite3_prepare_v2(db, insertFireshedsSQLString.c_str(), -1, &fireshedsStmt, NULL); // Prepare SQL statement
    if (rc != SQLITE_OK)
    {
        printf("Error: Could not create prepare SQL statement\n\n");
        return rc;
    }

    rc = sqlite3_prepare_v2(db, insertTotalsFiresSQLString.c_str(), -1, &totalsStmt, NULL); // Prepare SQL statement
    if (rc != SQLITE_OK)
    {
        printf("Error: Could not create prepare total fires SQL statement\n\n");
        return rc;
    }

    rc = sqlite3_exec(db, "PRAGMA SYNCHRONOUS=OFF", NULL, NULL, NULL);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &sqlErrMsg);

    if (verbose)
    {
        printf("Creating \"firesheds.db\", please wait...\n\n");
    }

    for (int wfipsCellIndex = 0; wfipsCellIndex < fireshedData.originCellsForWfipscell.size(); wfipsCellIndex++)
    {
        for (int originCellIndex = 0; originCellIndex < fireshedData.originCellsForWfipscell[wfipsCellIndex].size(); originCellIndex++)
        {
            wfipscell = fireshedData.finalIndexToWfipsCellMap.at(wfipsCellIndex);
            bindColumnIndex = sqlite3_bind_parameter_index(fireshedsStmt, ":wfipscell");
            rc = sqlite3_bind_int(fireshedsStmt, bindColumnIndex, wfipscell);

            origin = fireshedData.originCellsForWfipscell[wfipsCellIndex][originCellIndex];
            bindColumnIndex = sqlite3_bind_parameter_index(fireshedsStmt, ":origin");
            rc = sqlite3_bind_int(fireshedsStmt, bindColumnIndex, origin);

            numPairs = fireshedData.numWfipscellOriginPairs[wfipsCellIndex][originCellIndex];
            bindColumnIndex = sqlite3_bind_parameter_index(fireshedsStmt, ":num_pairs");
            rc = sqlite3_bind_int(fireshedsStmt, bindColumnIndex, numPairs);

            rc = sqlite3_step(fireshedsStmt);
            rc = sqlite3_reset(fireshedsStmt);
        }
    }

    fireshedData.finalIndexToWfipsCellMap.clear();
    fireshedData.originCellsForWfipscell.clear();
    fireshedData.numWfipscellOriginPairs.clear();

    rc = sqlite3_reset(fireshedsStmt);
    rc = sqlite3_finalize(fireshedsStmt);

    int totalForOrigin = -1;
    for (auto iterator = fireshedData.finalOriginCellToTotalPairCountMap.begin(); iterator != fireshedData.finalOriginCellToTotalPairCountMap.end(); iterator++)
    {
        origin = iterator->first;
        bindColumnIndex = sqlite3_bind_parameter_index(totalsStmt, ":origin");
        rc = sqlite3_bind_int(totalsStmt, bindColumnIndex, origin);

        totalForOrigin = iterator->second;
        bindColumnIndex = sqlite3_bind_parameter_index(totalsStmt, ":total_fires");
        rc = sqlite3_bind_int(totalsStmt, bindColumnIndex, totalForOrigin);

        rc = sqlite3_step(totalsStmt);
        rc = sqlite3_reset(totalsStmt);
    }

    rc = sqlite3_reset(totalsStmt);
    rc = sqlite3_finalize(totalsStmt);

    rc = sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, "PRAGMA SYNCHRONOUS=ON", NULL, NULL, NULL);

    //create indices
    printf("Creating indices on tables...\n");
    rc = sqlite3_exec(db, "CREATE INDEX idx_firesheds_wfipscell ON firesheds(wfipscell ASC)",
        NULL, NULL, &sqlErrMsg);

    rc = sqlite3_exec(db, "CREATE UNIQUE INDEX idx_totals_for_origins_origin ON totals_for_origins(origin ASC)",
        NULL, NULL, &sqlErrMsg);

    // "Vacuum" the database to free unused memory
    if (verbose)
    {
        printf("Optimizing database file size, please wait...\n\n");
    }
    rc = sqlite3_exec(db, "VACUUM", NULL, NULL, NULL);

    rc = sqlite3_close(db);

    sqlite3_free(sqlErrMsg);

    return rc;
}

bool AreClose(double a, double b)
{
    return fabs(a - b) < epsilon;
}

SBoundingBox GetBoundingBox(vector<MyPoint2D> theRingString)
{
    SBoundingBox theBoundingBox;
    theBoundingBox.minX = theRingString[0].X;
    theBoundingBox.maxX = theRingString[0].X;
    theBoundingBox.minY = theRingString[0].Y;
    theBoundingBox.maxY = theRingString[0].Y;

    for (int i = 0; i < theRingString.size(); i++)
    {
        if (theRingString[i].X < theBoundingBox.minX)
        {
            theBoundingBox.minX = theRingString[i].X;
        }
        if (theRingString[i].X > theBoundingBox.maxX)
        {
            theBoundingBox.maxX = theRingString[i].X;
        }
        if (theRingString[i].Y < theBoundingBox.minY)
        {
            theBoundingBox.minY = theRingString[i].Y;
        }
        if (theRingString[i].Y > theBoundingBox.maxY)
        {
            theBoundingBox.maxY = theRingString[i].Y;
        }
    }

    return theBoundingBox;
}

void Usage()
{
    printf("\nUsage:\n");
    printf("firesheds <--shape path>  Required\n");
    printf("          [--grid path] Optional\n");
    printf("          [--out path]  Optional\n");
    printf("          [--t number]  Optional\n");
    printf("          [--verbose]   Optional\n");
    printf("--shape <path>          Required: Specifies path to input shapefiles\n");
    printf("--grid <path>           Optional: Specifies path to WFIPSgrid.tif\n");
    printf("                        default: grid path is same as shape path\n");
    printf("--out <path>            Optional: Specifies output path\n");
    printf("                        default: output path is same as shape path\n");
    printf("--t <number>            Optional: Specifies number of threads\n");
    printf("                        default: number of available processors\n");
    printf("--verbose               Optional: Provides detailed progress information\n");
    printf("\n\nA valid path to FSIM shapefile data must exist\n");
    printf("\nA valid path to WFIPSgrid.tif must exist (same path as shapefiles by default)\n");
    printf("\nIf out path is specified, it must exist (same path as shapefiles by default)\n");
    printf("\nIf thread number is specified, it must be greater than zero\n");
    printf("Example usage:\n");
    printf("firesheds --shape C:/my_shapefiles --grid C:/grid --out C:/output --t 4 --verbose\n\n");

    exit(1); // Exit with error code 1
}
