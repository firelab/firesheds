#pragma once
#include <vector>
#include <string>
#include "gdal_priv.h"
#include "ogr_spatialref.h"
#include <Windows.h>

using namespace std;

typedef enum {WGT_UNKNOWN = 0, WGT_BOOL = 1, WGT_FLOAT = 2, WGT_INT = 3} WFIPS_GRID_TYPE;

class CWfipsGrid
{
public:
	CWfipsGrid();
	~CWfipsGrid();
	int LoadData(string fileName, WFIPS_GRID_TYPE gridType = WGT_UNKNOWN);
	bool Create(WFIPS_GRID_TYPE gridType, int numRows, int numCols, double res, double minX, double maxX, double minY, double maxY, double *pGeotransform, double ndValue);
	bool CellValueBool(int WfipsSRS, double lat, double lon);
	float CellValueFloat(int WfipsSRS, double lat, double lon);
	int CellValueInt(int WfipsSRS, double lat, double lon);
	double GetNoDataVal();
	bool CellValueBool(int row, int col);
	float CellValueFloat(int row, int col);
	int CellValueInt(int row, int col);

	int CellValueDirectInt(double lat, double lon);
	int CellValueDirectBool(double lat, double lon);
	int CellValueDirectFloat(double lat, double lon);

	bool SetCellValue(int WfipsSRS, double lat, double lon, bool val);
	bool SetCellValue(int WfipsSRS, double lat, double lon, int val);
	bool SetCellValue(int WfipsSRS, double lat, double lon, float val);
	bool SetCellValue(int row, int col, bool val);
	bool SetCellValue(int row, int col, int val);
	bool SetCellValue(int row, int col, float val);
	int GetNumX();
	int GetNumY();
	double GetRes();
	double GetNorth();
	double GetSouth();
	double GetEast();
	double GetWest();
	bool GetGeoTransform(double *pTargetTransform);
	bool SaveAs(string outfileName);
	OGRCoordinateTransformation *GetCoordinateTransformation(int SRS);
	OGRSpatialReference *GetSRS(int wfipsSRS);
private:
	long long int CellIndex(int WfipsSRS, double lat, double lon);
	//int GetWfipsGeotransform(double *pAdfGeoTransform);
	int SetGeoTransform(double *pAdfGeoTransform);
	string m_strFileName;
	WFIPS_GRID_TYPE m_gridType;
	bool *m_bVals;
	float *m_fVals;
	int *m_iVals;
	int m_nRows;
	int m_nCols;
	double m_minX;
	double m_maxX;
	double m_minY;
	double m_maxY;
	double m_cellResX;
	double m_cellResY;
	//gdal structures
	const char *m_pszWkt;
	double m_adfGeoTransform[6];
	double m_adfInvGeoTransform[6];
	OGRSpatialReference m_oSrcSRS;
	double m_dNodataVal;
	//coordinate transformations for the three WFIPS view settings
	OGRSpatialReference m_wfipsSRS[3];
	OGRCoordinateTransformation *m_wfipsCT[3];
	CRITICAL_SECTION TransformCS;
};
