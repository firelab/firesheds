#pragma once
#include "WfipsGrid.h"
#include <sqlite3.h>

class WfipsData;
//const __int64 NUM_EXCLUSION_ROWS = 29180;
//const __int64 NUM_EXCLUSION_COLS = 46260;
//const double EXCLUSION_RES = 100.0;
//const __int64 NUM_EXCLUSION_CELLS = NUM_EXCLUSION_ROWS * NUM_EXCLUSION_COLS;


//WFIPSGrid base layout
//X: 2313 Y: 1459 Bands: 1
//2000,-2000
//projection is EPSG:5070
//from WFIPSgrid.tif
// -2361582.2785416999831796,259071.7189354300498962 : 2264417.7214583000168204,3177071.7189354300498962
//-2361582.2785416999831796, 2264417.7214583000168204, 259071.7189354300498962, 3177071.7189354300498962

enum WG_LAYER{
	WGL_WFP = 0, WGL_AIR_TO_GROUND, WGL_DISC_SIZE, WGL_ESL_SIZE, WGL_ESL_TIME, WGL_FIRST_DELAY, WGL_EXCLUDED, WGL_PMP_ROLL_PERC,
	WGL_WALK_IN_BOAT, WGL_WALK_IN_CREW, WGL_WALK_IN_ENG, WGL_WALK_IN_HELITACK, WGL_WALK_IN_SMKJMP,
	WGL_WALK_IN_TRACKED, WGL_WALK_PERC, WGL_WATER_DROP, WGL_WILDERNESS, WGL_DISTROADS, WGL_EVACTIME, WGL_WUI, WGL_WATER, WGL_ELEC_GAS, WGL_END
};


class CWfipsGridData
{
	
public:
	CWfipsGridData();
	~CWfipsGridData();

	bool WG_GetCellCoords(int cell, double *minX, double *minY, double *maxX, double *maxY);
	const int WG_GetCellIndex(double x, double y);
	int GetWfipsCell(double lat, double lon, int WfipsSRS = 0);
	int GetWfipsCellDirect(double lat, double lon);
	int LoadData(const char* dataPath);
	int GetWfipsGridGeotransform(double adfGeotransform[6]);
	//int LoadExclusionData(const char* fileName);
	//bool GetExcluded(int WfipsSRS, double lat, double lon);
	CWfipsGrid *CreateNewWfipsGrid(WFIPS_GRID_TYPE gridType, double nodataVal);
	int AssignDisplocAreaWfipsCells(char *displocAreaName, sqlite3 *displocAreaDB);

	int GetNumX();
	int GetNumY();
	double GetRes();
	double GetNorth();
	double GetSouth();
	double GetEast();
	double GetWest();
	int GetWfipsCellIndex(int row, int col);
	int GetRowColFromWfipsCell(int cell, int *row, int *col);

	int GetGridValueAsInt(WG_LAYER layer, int WfipsSRS, double lat, double lon);
	float GetGridValueAsFloat(WG_LAYER layer, int WfipsSRS, double lat, double lon);
	bool GetGridValueAsBool(WG_LAYER layer, int WfipsSRS, double lat, double lon);
	double GetGridNoDataVal(WG_LAYER layer);

	bool GetExcluded(double lat, double lon, int WfipsSRS = 0);
	float GetAirToGround(double lat, double lon, int WfipsSRS = 0);
	float GetDiscSize(double lat, double lon, int WfipsSRS = 0);
	float GetDistRoads(double lat, double lon, int WfipsSRS = 0);
	int GetEslSize(double lat, double lon, int WfipsSRS = 0);
	int GetEslTime(double lat, double lon, int WfipsSRS = 0);
	int GetFirstDelay(double lat, double lon, int WfipsSRS = 0);
	int GetPumpRollPerc(double lat, double lon, int WfipsSRS = 0);
	int GetWalkInBoat(double lat, double lon, int WfipsSRS = 0);
	int GetWalkInCrew(double lat, double lon, int WfipsSRS = 0);
	int GetWalkInEng(double lat, double lon, int WfipsSRS = 0);
	int GetWalkInHelitack(double lat, double lon, int WfipsSRS = 0);
	int GetWalkInSmkJmp(double lat, double lon, int WfipsSRS = 0);
	int GetWalkInTracked(double lat, double lon, int WfipsSRS = 0);
	int GetWalkPerc(double lat, double lon, int WfipsSRS = 0);
	int GetWaterDrop(double lat, double lon, int WfipsSRS = 0);
	int GetElecGas(double lat, double lon, int WfipsSRS = 0);
	//added 2017/10/17 to accomodate calculating dispatch values
	//float GetDistRoads(double lat, double lon, int WfipsSRS = 0);
	int GetWater(double lat, double lon, int WfipsSRS = 0);
	int GetWilderness(double lat, double lon, int WfipsSRS = 0);
	int GetGroungEvacTimeClass(double lat, double lon, int WfipsSRS = 0);
	int GetWUIClass(double lat, double lon, int WfipsSRS = 0);

	bool GetExcluded(int row, int col);
	float GetAirToGround(int row, int col);
	float GetDiscSize(int row, int col);
	float GetDistRoads(int row, int col);
	int GetEslSize(int row, int col);
	int GetEslTime(int row, int col);
	int GetFirstDelay(int row, int col);
	int GetPumpRollPerc(int row, int col);
	int GetWalkInBoat(int row, int col);
	int GetWalkInCrew(int row, int col);
	int GetWalkInEng(int row, int col);
	int GetWalkInHelitack(int row, int col);
	int GetWalkInSmkJmp(int row, int col);
	int GetWalkInTracked(int row, int col);
	int GetWalkPerc(int row, int col);
	int GetWaterDrop(int row, int col);
	int GetElecGas(int row, int col);
	int GetWater(int row, int col);
	int GetWilderness(int row, int col);
	int GetGroungEvacTimeClass(int row, int col);
	int GetWUIClass(int row, int col);


	int GetDispatchValue(int dtIndex, double latitude, double longitude);
	int GetDispatchValue(int dtIndex, int wfipsGridRow, int wfipsGridCol);
	int GetDispatchValue(int dtIndex, double roadDist, int water, int WUI, int wild);
	char *GetLayerName(WG_LAYER layer);
	WFIPS_GRID_TYPE GetLayerType(WG_LAYER layer);
	char *GetLayerUnits(WG_LAYER layer);
	int GetLayerDecimals(WG_LAYER layer);
	int GetNumCellsMatching(WG_LAYER layer, bool findVal);
	int GetNumCellsMatching(WG_LAYER layer, float findVal);
	int GetNumCellsMatching(WG_LAYER layer, int findVal);
	vector<int> GetCellsMatching(WG_LAYER layer, bool findVal);
	//return 0 for not allowed, 1 for allowed, 2 for preferred
	vector<int> GetCellsMatching(WG_LAYER layer, float findVal);
	//return 0 for not allowed, 1 for allowed, 2 for preferred
	vector<int> GetCellsMatching(WG_LAYER layer, int findVal);
	CWfipsGrid *GetWfipsGrid();
	CWfipsGrid *GetFwaGrid();
private:
	string m_strDataPath;
	vector<CWfipsGrid *> m_WGgrids;
	//string m_strExclusionGridName;
	CWfipsGrid m_wgWfipscells;
	CWfipsGrid m_wgFwaGrid;
	//WFIPSGrid main attributes
	//-2361582.2785416999831796, 2264417.7214583000168204, 259071.7189354300498962, 3177071.7189354300498962
	const int m_WGnumCellsX = 2313;
	const int m_WGnumCellsY = 1459;
	const double m_WGnorth = 3177071.7189354300498962;
	const double m_WGsouth = 259071.7189354300498962;
	const double m_WGwest = -2361582.2785416999831796;
	const double m_WGeast = 2264417.7214583000168204;
	const double m_WGres = 2000.0;
	const char* m_WGproj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs";
	double m_WGadfGeoTransform[6];
};

