#include "WfipsGridData.h"

#include "gdal_priv.h"
#include "ogrsf_frmts.h"

char *WGL_LAYERNAMES[WGL_END] =
{
	"Wildfire Hazard Potential",
	"Air/Ground Ratio",
	"Discovery Size",
	"Escape Size",
	"Escape Time",
	"First Delay",
	"Exclusions",
	"Pump and Roll %",
	"Walk In Boat",
	"Walk In Crew",
	"Walk In Engine",
	"Walk In Helitack",
	"Walk In Smoke Jumper",
	"Walk In Tracked",
	"Walk In %",
	"Water Drops",
	"Wilderness",
	"Distance To Road",
	"Ground Evac Time",
	"WUI",
	"Water Presence",
	"Electic/Gas",
};

char *WGL_FILENAMES[WGL_END] =
{
	"Wildland_Fire_Potential",
	"air_to_ground_2km",
	"disc_size_2km",
	"esl_size_2km",
	"esl_time_2km",
	"first_delay_2km",
	"excluded_2km",
	"pmp_roll_perc_2km",
	"walk_in_boat_2km",
	"walk_in_crew_2km",
	"walk_in_eng_2km",
	"walk_in_helitack_2km",
	"walk_in_smkjmp_2km",
	"walk_in_tracked_2km",
	"walk_perc_2km",
	"water_drop_2km",
	"Wilderness_2km",//"WildCombined_2km",
	"DistRoad_2km",
	"GroundEvacTime_2km",
	"WUI_2km",
	"water_2km",
	"elecgas_2km",
};

WFIPS_GRID_TYPE WGL_WG_TYPES[WGL_END] =
{
	WGT_INT,
	WGT_FLOAT,
	WGT_FLOAT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_BOOL,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_FLOAT,
	WGT_FLOAT,
	WGT_INT,
	WGT_INT,
	WGT_INT,
	WGT_INT,

};

char *WGL_UNITS[WGL_END] =
{
	"Class",
	"Ratio",
	"Acres",
	"Acres",
	"Minutes",
	"Minutes",
	"Class",
	"Percent",
	"Minutes",
	"Minutes",
	"Minutes",
	"Minutes",
	"Minutes",
	"Minutes",
	"Percent",
	"Class",
	"Class",
	"Meters",
	"Class",
	"Class",
	"Class",
	"Class",
};

int WGL_DECIMALS[WGL_END] = 
{
	0,
	2,
	3,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	0,
};

CWfipsGridData::CWfipsGridData()
{
	//m_strExclusionGridName = "";
	//m_wgExclusion = NULL;
	//GDALDriver *poDriver = (GDALDriver *)GDALGetDriverByName("MEM");
	//GDALDataset *poWGdataset = poDriver->Create("", m_WGnumCellsX, m_WGnumCellsY, 1, GDT_CInt32, NULL);

}


CWfipsGridData::~CWfipsGridData()
{
	//if (m_wgExclusion != NULL)
	//	delete m_wgExclusion;
	if (m_WGgrids.size() > 0)
	{
		for (int g = 0; g < WGL_END; g++)
		{
			delete m_WGgrids[g];
			m_WGgrids[g] = NULL;
		}
		m_WGgrids.clear();
	}

}

int CWfipsGridData::GetNumX()
{
	return m_WGnumCellsX;
}

int CWfipsGridData::GetNumY()
{
	return m_WGnumCellsY;
}
double CWfipsGridData::GetRes()
{
	return m_WGres;
}
double CWfipsGridData::GetNorth()
{
	return m_WGnorth;
}
double CWfipsGridData::GetSouth()
{
	return m_WGsouth;
}
double CWfipsGridData::GetEast()
{
	return m_WGeast;
}
double CWfipsGridData::GetWest()
{
	return m_WGwest;
}

int CWfipsGridData::GetWfipsCellIndex(int row, int col)
{
	int ret = -1;
	if (row >= 0 && col >= 0 && row < m_WGnumCellsY && col < m_WGnumCellsX)
		ret = row * m_WGnumCellsX + col;
	return ret;
}

int CWfipsGridData::GetRowColFromWfipsCell(int cell, int *row, int *col)
{
	*row = cell / m_WGnumCellsX;
	*col = cell % m_WGnumCellsX;
	if (*row < 0)
		return -1;
	if(*row >= m_WGnumCellsY)
		return -2;
	if (*col < 0)
		return -3;
	if (*col >= m_WGnumCellsX)
		return -4;
	return 0;
}


int CWfipsGridData::LoadData(const char* dataPath)
{
	m_strDataPath = dataPath;
	string gridFileName;
	gridFileName = m_strDataPath + "WFIPSgrid.tif";
	m_wgWfipscells.LoadData(gridFileName, WGT_INT);
	m_wgWfipscells.GetGeoTransform(m_WGadfGeoTransform);
	//gridFileName = m_strDataPath + "\\fwa.tif";
	//m_wgFwaGrid.LoadData(gridFileName, WGT_INT);

	//for (int g = 0; g < WGL_END; g++)
	//{
	//	gridFileName = m_strDataPath + "\\" + WGL_FILENAMES[g] + ".tif";
	//	CWfipsGrid *pGrid = new CWfipsGrid();
	//	if(pGrid->LoadData(gridFileName, WGL_WG_TYPES[g]) == 0)
	//		m_WGgrids.push_back(pGrid);
	//	else
	//	{
	//		return -1;
	//	}
	//}
	return 0;
}

int CWfipsGridData::GetWfipsGridGeotransform(double adfGeotransform[6])
{
	for (int i = 0; i < 6; i++)
		adfGeotransform[i] = m_WGadfGeoTransform[i];
	return 0;
}

CWfipsGrid *CWfipsGridData::CreateNewWfipsGrid(WFIPS_GRID_TYPE gridType, double nodataVal)
{
	CWfipsGrid *pGrid = new CWfipsGrid();
	if (pGrid->Create(gridType, m_WGnumCellsY, m_WGnumCellsX, m_WGres, m_WGwest, m_WGeast, m_WGsouth, m_WGnorth, m_WGadfGeoTransform, nodataVal)
		!= true)
	{
		delete pGrid;
		pGrid = NULL;
	}
	return pGrid;
}

bool CWfipsGridData::WG_GetCellCoords(int cell, double *minX, double *minY, double *maxX, double *maxY)
{
	int col = cell % m_WGnumCellsX;
	int row = cell / m_WGnumCellsX;
	if (row < 0 || row >= m_WGnumCellsY
		|| col < 0 || col >= m_WGnumCellsX)
		return false;
	*minX = m_WGwest + col * m_WGres;
	*maxX = *minX + m_WGres;
	*maxY = m_WGnorth - row * m_WGres;
	*minY = *maxY - m_WGres;
	return true;
}

int CWfipsGridData::GetWfipsCell(double lat, double lon, int WfipsSRS/* = 0*/)
{
	return m_wgWfipscells.CellValueInt(WfipsSRS, lat, lon);
}

int CWfipsGridData::GetWfipsCellDirect(double lat, double lon)
{
	return m_wgWfipscells.CellValueDirectInt(lat, lon);
}

const int CWfipsGridData::WG_GetCellIndex(double x, double y)
{
	int ret = -1;
	if (x < m_WGwest || x > m_WGeast || y < m_WGsouth || y > m_WGnorth)
		return ret;
	//m_wgWfipscells.
	int col = (int)((x - m_WGwest) / m_WGres);
	int row = (int)((m_WGnorth - y) / m_WGres);
	return row * m_WGnumCellsX + col;
}

/*typedef struct
{
	int disploc_key;
	vector<int> cells;
} displocCellStruct;*/


int CWfipsGridData::AssignDisplocAreaWfipsCells(char *displocAreaName, sqlite3 *displocAreaDB)
{
	if (strlen(displocAreaName) <= 0)
		return 0;
	//CVectorLayer displocLayer;
	//displocLayer.Create("")

}

bool CWfipsGridData::GetExcluded(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_EXCLUDED])
	{
		return m_WGgrids[WGL_EXCLUDED]->CellValueBool(WfipsSRS, lat, lon);
	}
	return true;
}

float CWfipsGridData::GetAirToGround(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_AIR_TO_GROUND])
	{
		return m_WGgrids[WGL_AIR_TO_GROUND]->CellValueFloat(WfipsSRS, lat, lon);
	}
	return -9999.0;
}

float CWfipsGridData::GetDiscSize(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_DISC_SIZE])
	{
		return m_WGgrids[WGL_DISC_SIZE]->CellValueFloat(WfipsSRS, lat, lon);
	}
	return -9999.0;
}

int CWfipsGridData::GetEslSize(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_ESL_SIZE])
	{
		return m_WGgrids[WGL_ESL_SIZE]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetEslTime(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_ESL_TIME])
	{
		return m_WGgrids[WGL_ESL_TIME]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetFirstDelay(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_FIRST_DELAY])
	{
		return m_WGgrids[WGL_FIRST_DELAY]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetPumpRollPerc(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_PMP_ROLL_PERC])
	{
		return m_WGgrids[WGL_PMP_ROLL_PERC]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInBoat(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WALK_IN_BOAT])
	{
		return m_WGgrids[WGL_WALK_IN_BOAT]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInCrew(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WALK_IN_CREW])
	{
		return m_WGgrids[WGL_WALK_IN_CREW]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInEng(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WALK_IN_ENG])
	{
		return m_WGgrids[WGL_WALK_IN_ENG]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInHelitack(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WALK_IN_HELITACK])
	{
		return m_WGgrids[WGL_WALK_IN_HELITACK]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInSmkJmp(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WALK_IN_SMKJMP])
	{
		return m_WGgrids[WGL_WALK_IN_SMKJMP]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInTracked(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WALK_IN_TRACKED])
	{
		return m_WGgrids[WGL_WALK_IN_TRACKED]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWalkPerc(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WALK_PERC])
	{
		return m_WGgrids[WGL_WALK_PERC]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWaterDrop(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WATER_DROP])
	{
		return m_WGgrids[WGL_WATER_DROP]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetElecGas(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_ELEC_GAS])
	{
		return m_WGgrids[WGL_ELEC_GAS]->CellValueInt(WfipsSRS, lat, lon);
	}
	return true;
}

int CWfipsGridData::GetWater(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WATER])
	{
		return m_WGgrids[WGL_WATER]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

float CWfipsGridData::GetDistRoads(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_DISTROADS])
	{
		return m_WGgrids[WGL_DISTROADS]->CellValueFloat(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWilderness(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WILDERNESS])
	{
		return m_WGgrids[WGL_WILDERNESS]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetGroungEvacTimeClass(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_EVACTIME])
	{
		return m_WGgrids[WGL_EVACTIME]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

int CWfipsGridData::GetWUIClass(double lat, double lon, int WfipsSRS/* = 0*/)
{
	if (m_WGgrids[WGL_WUI])
	{
		return m_WGgrids[WGL_WUI]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;
}

bool CWfipsGridData::GetExcluded(int row, int col)
{
	if (m_WGgrids[WGL_EXCLUDED])
	{
		return m_WGgrids[WGL_EXCLUDED]->CellValueBool(row, col);
	}
	return true;
}

float CWfipsGridData::GetAirToGround(int row, int col)
{
	if (m_WGgrids[WGL_AIR_TO_GROUND])
	{
		return m_WGgrids[WGL_AIR_TO_GROUND]->CellValueFloat(row, col);
	}
	return -9999.0;
}

float CWfipsGridData::GetDiscSize(int row, int col)
{
	if (m_WGgrids[WGL_DISC_SIZE])
	{
		return m_WGgrids[WGL_DISC_SIZE]->CellValueFloat(row, col);
	}
	return -9999.0;
}

float CWfipsGridData::GetDistRoads(int row, int col)
{
	if (m_WGgrids[WGL_DISTROADS])
	{
		return m_WGgrids[WGL_DISTROADS]->CellValueFloat(row, col);
	}
	return -9999.0;
}

int CWfipsGridData::GetEslSize(int row, int col)
{
	if (m_WGgrids[WGL_ESL_SIZE])
	{
		return m_WGgrids[WGL_ESL_SIZE]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetEslTime(int row, int col)
{
	if (m_WGgrids[WGL_ESL_TIME])
	{
		return m_WGgrids[WGL_ESL_TIME]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetFirstDelay(int row, int col)
{
	if (m_WGgrids[WGL_FIRST_DELAY])
	{
		return m_WGgrids[WGL_FIRST_DELAY]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetPumpRollPerc(int row, int col)
{
	if (m_WGgrids[WGL_PMP_ROLL_PERC])
	{
		return m_WGgrids[WGL_PMP_ROLL_PERC]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInBoat(int row, int col)
{
	if (m_WGgrids[WGL_WALK_IN_BOAT])
	{
		return m_WGgrids[WGL_WALK_IN_BOAT]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInCrew(int row, int col)
{
	if (m_WGgrids[WGL_WALK_IN_CREW])
	{
		return m_WGgrids[WGL_WALK_IN_CREW]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInEng(int row, int col)
{
	if (m_WGgrids[WGL_WALK_IN_ENG])
	{
		return m_WGgrids[WGL_WALK_IN_ENG]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInHelitack(int row, int col)
{
	if (m_WGgrids[WGL_WALK_IN_HELITACK])
	{
		return m_WGgrids[WGL_WALK_IN_HELITACK]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInSmkJmp(int row, int col)
{
	if (m_WGgrids[WGL_WALK_IN_SMKJMP])
	{
		return m_WGgrids[WGL_WALK_IN_SMKJMP]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWalkInTracked(int row, int col)
{
	if (m_WGgrids[WGL_WALK_IN_TRACKED])
	{
		return m_WGgrids[WGL_WALK_IN_TRACKED]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWalkPerc(int row, int col)
{
	if (m_WGgrids[WGL_WALK_PERC])
	{
		return m_WGgrids[WGL_WALK_PERC]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWaterDrop(int row, int col)
{
	if (m_WGgrids[WGL_WATER_DROP])
	{
		return m_WGgrids[WGL_WATER_DROP]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetElecGas(int row, int col)
{
	if (m_WGgrids[WGL_ELEC_GAS])
	{
		return m_WGgrids[WGL_ELEC_GAS]->CellValueInt(row, col);
	}
	return true;
}

int CWfipsGridData::GetWater(int row, int col)
{
	if (m_WGgrids[WGL_WATER])
	{
		return m_WGgrids[WGL_WATER]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWilderness(int row, int col)
{
	if (m_WGgrids[WGL_WILDERNESS])
	{
		return m_WGgrids[WGL_WILDERNESS]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetGroungEvacTimeClass(int row, int col)
{
	if (m_WGgrids[WGL_EVACTIME])
	{
		return m_WGgrids[WGL_EVACTIME]->CellValueInt(row, col);
	}
	return -9999;
}

int CWfipsGridData::GetWUIClass(int row, int col)
{
	if (m_WGgrids[WGL_WUI])
	{
		return m_WGgrids[WGL_WUI]->CellValueInt(row, col);
	}
	return -9999;
}





int CWfipsGridData::GetGridValueAsInt(WG_LAYER layer, int WfipsSRS, double lat, double lon)
{
	if (m_WGgrids[layer])
	{
		return m_WGgrids[layer]->CellValueInt(WfipsSRS, lat, lon);
	}
	return -9999;

}

float CWfipsGridData::GetGridValueAsFloat(WG_LAYER layer, int WfipsSRS, double lat, double lon)
{
	if (m_WGgrids[layer])
	{
		return m_WGgrids[layer]->CellValueFloat(WfipsSRS, lat, lon);
	}
	return -9999;

}

bool CWfipsGridData::GetGridValueAsBool(WG_LAYER layer, int WfipsSRS, double lat, double lon)
{
	if (m_WGgrids[layer])
	{
		//switch (m_WGgrids[layer]->CellValueInt)
		return m_WGgrids[layer]->CellValueBool(WfipsSRS, lat, lon);
	}
	return -9999;

}

int CWfipsGridData::GetNumCellsMatching(WG_LAYER layer, bool findVal)
{
	int matches = 0;
	if (m_WGgrids[layer])
	{
		for (int r = 0; r < m_WGgrids[layer]->GetNumY(); r++)
		{
			for (int c = 0; c < m_WGgrids[layer]->GetNumX(); c++)
			{
				if (m_WGgrids[layer]->CellValueBool(r, c) == findVal)
					matches++;
			}
		}
	}
	return matches;
}

int CWfipsGridData::GetNumCellsMatching(WG_LAYER layer, float findVal)
{
	int matches = 0;
	if (m_WGgrids[layer])
	{
		for (int r = 0; r < m_WGgrids[layer]->GetNumY(); r++)
		{
			for (int c = 0; c < m_WGgrids[layer]->GetNumX(); c++)
			{
				if (m_WGgrids[layer]->CellValueFloat(r, c) == findVal)
					matches++;
			}
		}
	}
	return matches;
}

int CWfipsGridData::GetNumCellsMatching(WG_LAYER layer, int findVal)
{
	int matches = 0;
	if (m_WGgrids[layer])
	{
		for (int r = 0; r < m_WGgrids[layer]->GetNumY(); r++)
		{
			for (int c = 0; c < m_WGgrids[layer]->GetNumX(); c++)
			{
				if (m_WGgrids[layer]->CellValueInt(r, c) == findVal)
					matches++;
			}
		}
	}
	return matches;
}

vector<int> CWfipsGridData::GetCellsMatching(WG_LAYER layer, bool findVal)
{
	vector<int> cells;
	if (m_WGgrids[layer])
	{
		for (int r = 0; r < m_WGgrids[layer]->GetNumY(); r++)
		{
			for (int c = 0; c < m_WGgrids[layer]->GetNumX(); c++)
			{
				if (m_WGgrids[layer]->CellValueBool(r, c) == findVal)
					cells.push_back(GetWfipsCellIndex(r, c));
			}
		}
	}
	return cells;
}

vector<int> CWfipsGridData::GetCellsMatching(WG_LAYER layer, float findVal)
{
	vector<int> cells;
	if (m_WGgrids[layer])
	{
		for (int r = 0; r < m_WGgrids[layer]->GetNumY(); r++)
		{
			for (int c = 0; c < m_WGgrids[layer]->GetNumX(); c++)
			{
				if (m_WGgrids[layer]->CellValueFloat(r, c) == findVal)
					cells.push_back(GetWfipsCellIndex(r, c));
			}
		}
	}
	return cells;
}

vector<int> CWfipsGridData::GetCellsMatching(WG_LAYER layer, int findVal)
{
	vector<int> cells;
	if (m_WGgrids[layer])
	{
		for (int r = 0; r < m_WGgrids[layer]->GetNumY(); r++)
		{
			for (int c = 0; c < m_WGgrids[layer]->GetNumX(); c++)
			{
				if (m_WGgrids[layer]->CellValueInt(r, c) == findVal)
					cells.push_back(GetWfipsCellIndex(r, c));
			}
		}
	}
	return cells;
}

double CWfipsGridData::GetGridNoDataVal(WG_LAYER layer)
{
	return m_WGgrids[layer]->GetNoDataVal();
}

char *CWfipsGridData::GetLayerName(WG_LAYER layer)
{
	return WGL_LAYERNAMES[layer];
}

WFIPS_GRID_TYPE CWfipsGridData::GetLayerType(WG_LAYER layer)
{
	return WGL_WG_TYPES[layer];
}

char *CWfipsGridData::GetLayerUnits(WG_LAYER layer)
{
	return WGL_UNITS[layer];
}

int CWfipsGridData::GetLayerDecimals(WG_LAYER layer)
{
	return WGL_DECIMALS[layer];
}

//return 0 for not allowed, 1 for allowed, 2 for preferred
int CWfipsGridData::GetDispatchValue(int dtIndex, double latitude, double longitude)
{
	//sets the dispatch levels for the fire based on WfipsGridData
	//vareis by class, but mostly uses DistRoads() layer
	//secondary processing on boats based on presence of water
	//secondary processing on ground equipment (Dozers, Engines, etc) based on Wilderness/Roadless
	//uses global wfipsGridData to get cell values based on lat/lon of fire
	//levels are:
	//	0 = not allowed
	//	1 = allowed
	//	2 = preferred
	double roadDist = GetDistRoads(latitude, longitude);
	int water = GetWater(latitude, longitude);
	int WUI = GetWUIClass(latitude, longitude);
	int wild = GetWilderness(latitude, longitude);
	return GetDispatchValue(dtIndex, roadDist, water, WUI, wild);

}

//return 0 for not allowed, 1 for allowed, 2 for preferred
int CWfipsGridData::GetDispatchValue(int dtIndex, int wfipsGridRow, int wfipsGridCol)
{
	//levels are:
	//	0 = not allowed
	//	1 = allowed
	//	2 = preferred
	double roadDist = GetDistRoads(wfipsGridRow, wfipsGridCol);
	int water = GetWater(wfipsGridRow, wfipsGridCol);
	int WUI = GetWUIClass(wfipsGridRow, wfipsGridCol);
	int wild = GetWilderness(wfipsGridRow, wfipsGridCol);
	return GetDispatchValue(dtIndex, roadDist, water, WUI, wild);
}

int CWfipsGridData::GetDispatchValue(int dtIndex, double roadDist, int water, int WUI, int wild)
{
	//levels are:
	//	0 = not allowed
	//	1 = allowed
	//	2 = preferred
	int val = 0;
	//switch (dtIndex)
	//{
	//case DTI_ATT://airtankers
	//	val = 1;
	//	if (WUI > 1)
	//		val = 2;
	//	else if (wild > 0)
	//		val = 0;
	//	break;
	//case DTI_CRW://crews
	//	val = 1;
	//	if (roadDist < 2000.0)//2 km
	//		val = 2; break;
	//case DTI_DZR://dozers
	//	val = 1;
	//	if (roadDist < 200.0 && roadDist >= 0) //clost to road, prefered
	//		val = 2;
	//	if (wild > 0) //no dozers in wilderness
	//		val = 0;
	//	if (roadDist > 1000)// > 1Km no dozers
	//		val = 0;
	//	break;
	//case DTI_EN://engines
	//	val = 1;
	//	if (roadDist < 500.0 && roadDist >= 0) //clost to road, prefered
	//		val = 2;
	//	break;
	//case DTI_FBT://boats
	//	val = 0;
	//	if (water > 0)
	//		val = 1;
	//	break;
	//case DTI_HEL:
	//	val = 1;
	//	if (roadDist >= 2000.0)
	//		val = 2;
	//	break;
	//case DTI_HELI:
	//	val = 1;
	//	if (roadDist >= 2000.0)
	//		val = 2;
	//	break;
	//case DTI_SEAT:
	//	val = 1;
	//	if (WUI > 1)
	//		val = 2;
	//	else if (wild > 0)
	//		val = 0;
	//	break;
	//case DTI_SCP:
	//	val = 1;
	//	if (water > 0)
	//		val = 2;
	//	if (wild > 0)
	//		val = 0;
	//	break;
	//case DTI_SMJR:
	//	val = 1;
	//	if (roadDist >= 2000.0)
	//		val = 2;
	//	break;
	//case DTI_TP:
	//	val = 1;
	//	if (roadDist < 200.0 && roadDist >= 0) //clost to road, prefered
	//		val = 2;
	//	if (wild > 0) //no TP in wilderness
	//		val = 0;
	//	if (roadDist > 1000)// > 1Km no TP
	//		val = 0;
	//	break;
	//case DTI_WT:
	//	val = 1;
	//	if (roadDist > 200.0)//200 meters don't use
	//		val = 0;
	//	break;
	//}
	return val;
}

CWfipsGrid *CWfipsGridData::GetWfipsGrid()
{
	return &m_wgWfipscells;
}

CWfipsGrid *CWfipsGridData::GetFwaGrid()
{
	return &m_wgFwaGrid;
}
