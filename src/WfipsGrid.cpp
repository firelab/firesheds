#include "WfipsGrid.h"

char *MapESRIProjStrings[] =
{
	"",
	"PROJCS[\"Albers\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Albers\"],PARAMETER[\"standard_parallel_1\",29.5],PARAMETER[\"standard_parallel_2\",45.5],PARAMETER[\"latitude_of_origin\",23],PARAMETER[\"central_meridian\",-96],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"METERS\",1]]",
	"PROJCS[\"WGS 84 / Pseudo - Mercator\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Mercator\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]"
};

CWfipsGrid::CWfipsGrid()
{
	m_pszWkt = "";
	m_strFileName = "";
	m_bVals = NULL;
	m_fVals = NULL;
	m_iVals = NULL;
	m_wfipsCT[0] = m_wfipsCT[1] = m_wfipsCT[2] = NULL;
	//InitializeCriticalSection(&TransformCS);
}


CWfipsGrid::~CWfipsGrid()
{
	if (m_bVals)
		delete[] m_bVals;
	if (m_fVals)
		delete[] m_fVals;
	if (m_iVals)
		delete[] m_iVals;
	if (m_wfipsCT[0])
		OGRCoordinateTransformation::DestroyCT(m_wfipsCT[0]);
	if (m_wfipsCT[1])
		OGRCoordinateTransformation::DestroyCT(m_wfipsCT[1]);
	if (m_wfipsCT[2])
		OGRCoordinateTransformation::DestroyCT(m_wfipsCT[2]);
	//DeleteCriticalSection(&TransformCS);

}

int CWfipsGrid::LoadData(string fileName, WFIPS_GRID_TYPE gridType)
{
	m_strFileName = fileName;
	GDALDataset  *poDataset;
	poDataset = (GDALDataset *)GDALOpen(m_strFileName.c_str(), GA_ReadOnly);
	if (poDataset == NULL)
	{
		m_strFileName = "";
		return -1;
	}

	m_pszWkt = poDataset->GetProjectionRef();
	int err = poDataset->GetGeoTransform(m_adfGeoTransform);
	if (err != CE_None)
	{
		GDALClose(poDataset);
		return -1;
	}
	err = GDALInvGeoTransform(m_adfGeoTransform, m_adfInvGeoTransform);
	err = m_oSrcSRS.importFromWkt((char**)&m_pszWkt);
	err = m_oSrcSRS.morphFromESRI();
	m_nRows = poDataset->GetRasterYSize();
	m_nCols = poDataset->GetRasterXSize();
	m_cellResX = fabs(m_adfGeoTransform[1]);
	m_cellResY = fabs(m_adfGeoTransform[5]);
	if (m_cellResX <= 0.0)
		m_cellResX = 1.0;
	if (m_cellResY <= 0.0)
		m_cellResY = 1.0;
	m_minX = m_adfGeoTransform[0];
	m_maxX = m_minX + m_nCols * m_cellResX;
	m_maxY = m_adfGeoTransform[3];
	m_minY = m_maxY - m_nRows * m_cellResY;
	GDALRasterBand *pBand = NULL;
	pBand = poDataset->GetRasterBand(1);
	GDALDataType geoType = pBand->GetRasterDataType();
	m_dNodataVal = pBand->GetNoDataValue();
	void *panScanline = NULL;
	if (gridType != WGT_UNKNOWN)

		m_gridType = gridType;
	else
	{
		switch (geoType)
		{
		case GDT_Byte:
		case GDT_UInt16:
		case GDT_Int16:
		case GDT_UInt32:
		case GDT_Int32:
			m_gridType = WGT_INT;
			break;
		case GDT_Float32:
		case GDT_Float64:
			m_gridType = WGT_FLOAT;
			break;
		default:
			m_gridType = WGT_INT;
			break;

		}
	}
	switch (m_gridType)
	{
	case WGT_BOOL:
		m_bVals = new bool[m_nRows * m_nCols];
		break;
	case WGT_INT:
		m_iVals = new int[m_nRows * m_nCols];
		break;
	case WGT_FLOAT:
		m_fVals = new float[m_nRows * m_nCols];
		break;
	}
	//fetch the data into proper vector
	switch (geoType)
	{
	case GDT_Byte:
		panScanline = (char*)CPLMalloc(sizeof(char) * m_nCols);
		break;
	case GDT_UInt16:
		panScanline = (unsigned short*)CPLMalloc(sizeof(unsigned short) * m_nCols);
		break;
	case GDT_Int16:
		panScanline = (short*)CPLMalloc(sizeof(short) * m_nCols);
		break;
	case GDT_UInt32:
		panScanline = (unsigned int*)CPLMalloc(sizeof(unsigned int) * m_nCols);
		break;
	case GDT_Int32:
		panScanline = (int*)CPLMalloc(sizeof(int) * m_nCols);
		break;
	case GDT_Float32:
		panScanline = (float*)CPLMalloc(sizeof(float) * m_nCols);
		break;
	case GDT_Float64:
		panScanline = (long long int*)CPLMalloc(sizeof(long long int) * m_nCols);
		break;
	default:
		break;
	}

	for (int i = 0; i < m_nRows; i++)
	{
		pBand->RasterIO(GF_Read, 0, i, m_nCols, 1, panScanline, m_nCols,
			1, geoType, 0, 0);
		for (int x = 0; x < m_nCols; x++)
		{
			switch (m_gridType)
			{
			case WGT_BOOL:
				switch (geoType)
				{
				case GDT_Byte:
					m_bVals[i * m_nCols + x] = (bool)((char *)panScanline)[x];
					break;
				case GDT_UInt16:
					m_bVals[i * m_nCols + x] = (bool)((unsigned short *)panScanline)[x];
					break;
				case GDT_Int16:
					m_bVals[i * m_nCols + x] = (bool)((short *)panScanline)[x];
					break;
				case GDT_UInt32:
					m_bVals[i * m_nCols + x] = (bool)((unsigned int *)panScanline)[x];
					break;
				case GDT_Int32:
					m_bVals[i * m_nCols + x] = (bool)((int *)panScanline)[x];
					break;
				case GDT_Float32:
					m_bVals[i * m_nCols + x] = (bool)((float *)panScanline)[x];
					break;
				case GDT_Float64:
					m_bVals[i * m_nCols + x] = (bool)((double *)panScanline)[x];
					break;
				default:
					m_bVals[i * m_nCols + x] = false;
					break;
				}
				break;
			case WGT_FLOAT:
				switch (geoType)
				{
				case GDT_Byte:
					m_fVals[i * m_nCols + x] = (float)((char *)panScanline)[x];
					break;
				case GDT_UInt16:
					m_fVals[i * m_nCols + x] = (float)((unsigned short *)panScanline)[x];
					break;
				case GDT_Int16:
					m_fVals[i * m_nCols + x] = (float)((short *)panScanline)[x];
					break;
				case GDT_UInt32:
					m_fVals[i * m_nCols + x] = (float)((unsigned int *)panScanline)[x];
					break;
				case GDT_Int32:
					m_fVals[i * m_nRows + x] = (float)((int *)panScanline)[x];
					break;
				case GDT_Float32:
					m_fVals[i * m_nCols + x] = (float)((float *)panScanline)[x];
					break;
				case GDT_Float64:
					m_fVals[i * m_nCols + x] = (float)((double *)panScanline)[x];
					break;
				default:
					m_fVals[i * m_nCols + x] = m_dNodataVal;
					break;
				}
				break;
			case WGT_INT:
				switch (geoType)
				{
				case GDT_Byte:
					m_iVals[i * m_nCols + x] = (int)((char *)panScanline)[x];
					break;
				case GDT_UInt16:
					m_iVals[i * m_nCols + x] = (int)((unsigned short *)panScanline)[x];
					break;
				case GDT_Int16:
					m_iVals[i * m_nCols + x] = (int)((short *)panScanline)[x];
					break;
				case GDT_UInt32:
					m_iVals[i * m_nCols + x] = (int)((unsigned int *)panScanline)[x];
					break;
				case GDT_Int32:
					m_iVals[i * m_nCols + x] = (int)((int *)panScanline)[x];
					break;
				case GDT_Float32:
					m_iVals[i * m_nCols + x] = (int)((float *)panScanline)[x];
					break;
				case GDT_Float64:
					m_iVals[i * m_nCols + x] = (int)((double *)panScanline)[x];
					break;
				default:
					m_iVals[i * m_nCols + x] = (int)m_dNodataVal;
					break;
				}
				break;
			}
		}
	}
	CPLFree(panScanline);
	err = m_wfipsSRS[0].SetWellKnownGeogCS("NAD83");
	CPLStringList papszPrj5070 = NULL;
	papszPrj5070.AddString(MapESRIProjStrings[1]);
	err = m_wfipsSRS[1].importFromESRI(papszPrj5070);
	//err = m_wfipsSRS[1].SetWellKnownGeogCS("EPSG:5070");
	CPLStringList papszPrj3857 = NULL;
	papszPrj3857.AddString(MapESRIProjStrings[2]);
	err = m_wfipsSRS[2].importFromESRI(papszPrj3857);

	//err = m_wfipsSRS[2].SetWellKnownGeogCS("EPSG:3857");
	m_wfipsCT[0] = OGRCreateCoordinateTransformation(&m_wfipsSRS[0], &m_oSrcSRS);
	m_wfipsCT[1] = OGRCreateCoordinateTransformation(&m_wfipsSRS[1], &m_oSrcSRS);
	m_wfipsCT[2] = OGRCreateCoordinateTransformation(&m_wfipsSRS[2], &m_oSrcSRS);
	GDALClose(poDataset);
	return 0;
}

bool CWfipsGrid::SaveAs(string outfileName)
{
	GDALDriver *poTiffDriver;
	GDALDataset *poDstDS;
	/*
	* Apply DEFLATE compression to the file and use a geotiff profile
	*/
	char **papszCreateOptions = NULL;
	papszCreateOptions = CSLAddString(papszCreateOptions, "COMPRESS=DEFLATE");
	papszCreateOptions = CSLAddString(papszCreateOptions, "PROFILE=GDALGeoTIFF");

	poTiffDriver = (GDALDriver*)GDALGetDriverByName("GTiff");

	GDALDataType outType = GDT_Int32;
	if (m_gridType == WGT_BOOL)
		outType = GDT_Byte;
	else if (m_gridType == WGT_FLOAT)
		outType = GDT_Float32;

	poDstDS = (GDALDataset*)poTiffDriver->Create(outfileName.c_str(), m_nCols, m_nRows,
		1, outType,
		papszCreateOptions);
	if (poDstDS == NULL)
	{
		m_strFileName = "";
		return false;
	}

	poDstDS->SetProjection(MapESRIProjStrings[1]);
	poDstDS->SetGeoTransform(m_adfGeoTransform);
	poDstDS->SetDescription("WFIPS output GeoTIFF");

	GDALRasterBand *poDstBand;
	poDstBand = poDstDS->GetRasterBand(1);
	poDstBand->SetDescription("Triggered WFIPScells");
	int err;
	void *panScanLine;
	switch (outType)
	{
	case GDT_Int32:
		panScanLine = (int*)CPLMalloc(sizeof(int) * m_nCols);
		break;
	case GDT_Byte:
		panScanLine = (char*)CPLMalloc(sizeof(char) * m_nCols);
		break;
	case GDT_Float32:
		panScanLine = (float*)CPLMalloc(sizeof(float) * m_nCols);
		break;
	}
	for (int r = 0; r < m_nRows; r++)
	{
		for (int c = 0; c < m_nCols; c++)
		{
			switch (outType)
			{
			case GDT_Byte:
				((char *)panScanLine)[c] = (char) m_bVals[r * m_nCols + c];
				break;
			case GDT_Int32:
				((int *)panScanLine)[c] = m_iVals[r * m_nCols + c];
				break;
			case GDT_Float32:
				((float *)panScanLine)[c] = m_fVals[r * m_nCols + c];
				break;
			}
		}
		err = poDstBand->RasterIO(GF_Write, 0, r, m_nCols, 1, panScanLine, m_nCols,
			1, outType, 0, 0);

	}
	err = poDstBand->SetNoDataValue(m_dNodataVal);

	CPLFree(panScanLine);
	CSLDestroy(papszCreateOptions);
	GDALClose(poDstDS);
	return true;
}

bool CWfipsGrid::Create(WFIPS_GRID_TYPE gridType, int numRows, int numCols, double res, double minX, double maxX, double minY, double maxY, double *pGeotransform, double ndValue)
{
	m_strFileName = "";
	m_gridType = gridType;
	m_nRows = numRows;
	m_nCols = numCols;
	m_cellResX = m_cellResY = res;
	m_minX = minX;
	m_maxX = maxX;
	m_minY = minY;
	m_maxY = maxY;
	switch (gridType)
	{
	case WGT_BOOL:
		m_bVals = new bool[m_nRows * m_nCols];
		memset(m_bVals, (bool)ndValue, m_nRows * m_nCols * sizeof(bool));
		break;
	case WGT_INT:
		m_iVals = new int[m_nRows * m_nCols];
		memset(m_iVals, (int)ndValue, m_nRows * m_nCols * sizeof(int));
		break;
	case WGT_FLOAT:
		m_fVals = new float[m_nRows * m_nCols];
		memset(m_fVals, (float)ndValue, m_nRows * m_nCols * sizeof(float));
		break;
	}
	m_dNodataVal = ndValue;

	int err = m_wfipsSRS[0].SetWellKnownGeogCS("NAD83");
	CPLStringList papszPrj5070 = NULL;
	papszPrj5070.AddString(MapESRIProjStrings[1]);
	err = m_wfipsSRS[1].importFromESRI(papszPrj5070);
	//err = m_wfipsSRS[1].SetWellKnownGeogCS("EPSG:5070");
	CPLStringList papszPrj3857 = NULL;
	papszPrj3857.AddString(MapESRIProjStrings[2]);
	err = m_wfipsSRS[2].importFromESRI(papszPrj3857);
	//err = m_wfipsSRS[2].SetWellKnownGeogCS("EPSG:3857");
	papszPrj5070 = NULL;
	papszPrj5070.AddString(MapESRIProjStrings[1]);
	err = m_oSrcSRS.importFromESRI(papszPrj5070);
	//coordinate transformations
	m_wfipsCT[0] = OGRCreateCoordinateTransformation(&m_wfipsSRS[0], &m_oSrcSRS);
	m_wfipsCT[1] = OGRCreateCoordinateTransformation(&m_wfipsSRS[1], &m_oSrcSRS);
	m_wfipsCT[2] = OGRCreateCoordinateTransformation(&m_wfipsSRS[2], &m_oSrcSRS);

	SetGeoTransform(pGeotransform);

	return true;
}

int CWfipsGrid::SetGeoTransform(double *pAdfGeoTransform)
{
	for (int i = 0; i < 6; i++)
		m_adfGeoTransform[i] = pAdfGeoTransform[i];
	int err = GDALInvGeoTransform(m_adfGeoTransform, m_adfInvGeoTransform);
	return err;
}

double CWfipsGrid::GetNoDataVal()
{
	return m_dNodataVal;
}

bool CWfipsGrid::CellValueBool(int row, int col)
{
    long long int loc = row * m_nCols + col;
	if (loc >= 0 && loc < m_nCols * m_nRows && m_bVals != NULL)
		return m_bVals[loc];
	return false;
}

float CWfipsGrid::CellValueFloat(int row, int col)
{
    long long int loc = row * m_nCols + col;
	if (loc >= 0 && loc < m_nCols * m_nRows && m_fVals != NULL)
		return m_fVals[loc];
	return m_dNodataVal;
}

int CWfipsGrid::CellValueInt(int row, int col)
{
    long long int loc = row * m_nCols + col;
	if (loc >= 0 && loc < m_nCols * m_nRows && m_iVals != NULL)
		return m_iVals[loc];
	return m_dNodataVal;
}

int CWfipsGrid::CellValueDirectInt(double lat, double lon)
{
	int nCol = (int)(m_adfInvGeoTransform[0] + m_adfInvGeoTransform[1] *
		lon + m_adfInvGeoTransform[2] * lat);
	int nRow = (int)(m_adfInvGeoTransform[3] + m_adfInvGeoTransform[4] *
		lon + m_adfInvGeoTransform[5] * lat);
	if (nRow >= 0 && nRow < m_nRows
		&& nCol >= 0 && nCol < m_nCols)
	{
		switch (m_gridType)
		{
		case WGT_BOOL:
			return (int)m_bVals[nRow * m_nCols + nCol];
			break;
		case WGT_INT:
			return m_iVals[nRow * m_nCols + nCol];
			break;
		case WGT_FLOAT:
			return (int)m_fVals[nRow * m_nCols + nCol];
			break;
		}
	}
	return m_dNodataVal;
}

int CWfipsGrid::CellValueDirectBool(double lat, double lon)
{
	int nCol = (int)(m_adfInvGeoTransform[0] + m_adfInvGeoTransform[1] *
		lon + m_adfInvGeoTransform[2] * lat);
	int nRow = (int)(m_adfInvGeoTransform[3] + m_adfInvGeoTransform[4] *
		lon + m_adfInvGeoTransform[5] * lat);
	if (nRow >= 0 && nRow < m_nRows
		&& nCol >= 0 && nCol < m_nCols)
	{
		switch (m_gridType)
		{
		case WGT_BOOL:
			return m_bVals[nRow * m_nCols + nCol];
			break;
		case WGT_INT:
			return (bool)m_iVals[nRow * m_nCols + nCol];
			break;
		case WGT_FLOAT:
			return (bool)m_fVals[nRow * m_nCols + nCol];
			break;
		}
	}
	return (bool)m_dNodataVal;
}

int CWfipsGrid::CellValueDirectFloat(double lat, double lon)
{
	int nCol = (int)(m_adfInvGeoTransform[0] + m_adfInvGeoTransform[1] *
		lon + m_adfInvGeoTransform[2] * lat);
	int nRow = (int)(m_adfInvGeoTransform[3] + m_adfInvGeoTransform[4] *
		lon + m_adfInvGeoTransform[5] * lat);
	if (nRow >= 0 && nRow < m_nRows
		&& nCol >= 0 && nCol < m_nCols)
	{
		switch (m_gridType)
		{
		case WGT_BOOL:
			return (float)m_bVals[nRow * m_nCols + nCol];
			break;
		case WGT_INT:
			return (float)m_iVals[nRow * m_nCols + nCol];
			break;
		case WGT_FLOAT:
			return m_fVals[nRow * m_nCols + nCol];
			break;
		}
	}
	return (float)m_dNodataVal;
}

bool CWfipsGrid::SetCellValue(int WfipsSRS, double lat, double lon, bool val)
{
    long long int loc = CellIndex(WfipsSRS, lat, lon);
	if (loc >= 0 && loc < m_nRows * m_nCols)
	{//	set appropriate data array value
		switch (m_gridType)
		{
		case WGT_BOOL:
			m_bVals[loc] = val;
		case WGT_INT:
			m_iVals[loc] = (int)val;
		case WGT_FLOAT:
			m_fVals[loc] = (float)val;
		}
		return true;
	}
	return false;
}

bool CWfipsGrid::SetCellValue(int WfipsSRS, double lat, double lon, int val)
{
    long long int loc = CellIndex(WfipsSRS, lat, lon);
	if (loc >= 0 && loc < m_nRows * m_nCols)
	{//	set appropriate data array value
		switch (m_gridType)
		{
		case WGT_BOOL:
			m_bVals[loc] = (bool)val;
		case WGT_INT:
			m_iVals[loc] = val;
		case WGT_FLOAT:
			m_fVals[loc] = (float)val;
		}
		return true;
	}
	return false;
}

bool CWfipsGrid::SetCellValue(int WfipsSRS, double lat, double lon, float val)
{
    long long int loc = CellIndex(WfipsSRS, lat, lon);
	if (loc >= 0 && loc < m_nRows * m_nCols)
	{//	set appropriate data array value
		switch (m_gridType)
		{
		case WGT_BOOL:
			m_bVals[loc] = (bool)val;
		case WGT_INT:
			m_iVals[loc] = (int)val;
		case WGT_FLOAT:
			m_fVals[loc] = val;
		}
		return true;
	}
	return false;
}

bool CWfipsGrid::SetCellValue(int row, int col, bool val)
{
    long long int loc = row * m_nCols + col;
	if (loc >= 0 && loc < m_nCols * m_nRows && m_bVals != NULL)
		m_bVals[loc] = val;
	return false;
}

bool CWfipsGrid::SetCellValue(int row, int col, int val)
{
    long long int loc = row * m_nCols + col;
	if (loc >= 0 && loc < m_nCols * m_nRows && m_iVals != NULL)
		m_iVals[loc] = val;
	return false;
}

bool CWfipsGrid::SetCellValue(int row, int col, float val)
{
    long long int loc = row * m_nCols + col;
	if (loc >= 0 && loc < m_nCols * m_nRows && m_fVals != NULL)
		m_fVals[loc] = val;
	return false;
}

int CWfipsGrid::GetNumX()
{
	return m_nCols;
}

int CWfipsGrid::GetNumY()
{
	return m_nRows;
}

double CWfipsGrid::GetRes()
{
	return m_cellResX;
}

double CWfipsGrid::GetNorth()
{
	return m_maxY;
}

double CWfipsGrid::GetSouth()
{
	return m_minY;
}

double CWfipsGrid::GetEast()
{
	return m_maxX;
}

double CWfipsGrid::GetWest()
{
	return m_minX;
}

bool CWfipsGrid::GetGeoTransform(double *pTargetTransform)
{
	for (int i = 0; i < 6; i++)
	{
		pTargetTransform[i] = m_adfGeoTransform[i];
	}
	return true;
}

//CellIndex for lat/long assumed to be wfips native NAD83 (EPSG:2639)
long long int CWfipsGrid::CellIndex(int WfipsSRS, double lat, double lon)
{
	if (m_wfipsCT[WfipsSRS] == NULL)
	{
		return -3;
	}
	double x = lon, y = lat;
	//EnterCriticalSection(&TransformCS);
	m_wfipsCT[WfipsSRS]->Transform(1, &x, &y);
	//LeaveCriticalSection(&TransformCS);
	int nCol = (int)(m_adfInvGeoTransform[0] + m_adfInvGeoTransform[1] *
		x + m_adfInvGeoTransform[2] * y);
	int nRow = (int)(m_adfInvGeoTransform[3] + m_adfInvGeoTransform[4] *
		x + m_adfInvGeoTransform[5] * y);
	if (nRow >= 0 && nRow < m_nRows
		&& nCol >= 0 && nCol < m_nCols)
		return nRow * m_nCols + nCol;
	return -1;
}

bool CWfipsGrid::CellValueBool(int WfipsSRS, double lat, double lon)
{
    long long int loc = CellIndex(WfipsSRS, lat, lon);
	if (loc >= 0 && loc < m_nRows * m_nCols)
	{//	return m_bVals[loc];
		switch (m_gridType)
		{
		case WGT_BOOL:
			return m_bVals[loc];
		case WGT_INT:
			return (bool)m_iVals[loc];
		case WGT_FLOAT:
			return (bool)m_fVals[loc];
		}
	}
	return false;
}

float CWfipsGrid::CellValueFloat(int WfipsSRS, double lat, double lon)
{
    long long int loc = CellIndex(WfipsSRS, lat, lon);
	if (loc >= 0 && loc < m_nRows * m_nCols)
	{//	return m_bVals[loc];
		switch (m_gridType)
		{
		case WGT_BOOL:
			return (float) m_bVals[loc];
		case WGT_INT:
			return (float)m_iVals[loc];
		case WGT_FLOAT:
			return m_fVals[loc];
		}
	}
	return m_dNodataVal;
}

int CWfipsGrid::CellValueInt(int WfipsSRS, double lat, double lon)
{
    long long int loc = CellIndex(WfipsSRS, lat, lon);
	if (loc >= 0 && loc < m_nRows * m_nCols)
	{
		switch (m_gridType)
		{
		case WGT_BOOL:
			return (int) m_bVals[loc];
		case WGT_INT:
			return m_iVals[loc];
		case WGT_FLOAT:
			return (int)m_fVals[loc];
		}
	}
	return m_dNodataVal;
}

OGRCoordinateTransformation *CWfipsGrid::GetCoordinateTransformation(int SRS)
{
	return m_wfipsCT[SRS];
}
