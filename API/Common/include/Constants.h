//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Constants.h
//! @brief     define some global const value.
//! @author    He Xin.

#pragma once
#include "Precision.h"

namespace PHSPACE
{
const int SOLVER_BASED        = 0;
const int GRID_BASED          = 1;

const int SOLVE_FIELD         = 0;
const int CREATE_GRID         = 1;
const int CAL_WALL_DIST       = 2;
const int PARTITION_GRID      = 3;
const int REPOSITORY          = 4;
const int OVERSETGRID_VIEW    = 5;
const int OVERSET_CONFIG      = 6;
const int CONVERT_PARTITION   = 7;
const int DSMC_SOLVER         = 8;
const int HELICOPTER          = 9;
const int UGKS_SOLVER         = 10;
const int STRPATCHGRID_CONFIG = 11;
const int HO_SOLVER           = 12;
const int LBM_SOLVER          = 13;

const int INTEGRATIVE_SOLVER  = 14;

const int SPECDIFFHYB_Solver = 15;

const int POST_PROCESSING     = 99;

const RDouble half            = 0.5;
const RDouble third           = 1.0 / 3.0;
const RDouble two3rd          = 2.0 / 3.0;
const RDouble four3rd         = 4.0 / 3.0;
const RDouble fourth          = 0.25;
const RDouble sixth           = 1.0 / 6.0;
const RDouble eighth          = 1.0 / 8.0;
const RDouble tenth           = 0.1;
const RDouble zero            = 0.0;
const RDouble one             = 1.0;
const RDouble two             = 2.0;
const RDouble three           = 3.0;
const RDouble four            = 4.0;
const RDouble five            = 5.0;
const RDouble six             = 6.0;
const RDouble seven           = 7.0;
const RDouble eight           = 8.0;
const RDouble nine            = 9.0;
const RDouble ten             = 10.0;
const RDouble eleven          = 11.0;
const RDouble twelve          = 12.0;
const RDouble thirteen        = 13.0;
const RDouble fourteen        = 14.0;
const RDouble fifteen         = 15.0;
const RDouble sixteen         = 16.0;
const RDouble seventeen       = 17.0;
const RDouble eighteen        = 18.0;
const RDouble nineteen        = 19.0;
const RDouble twenty          = 20.0;

const RDouble PI              = 3.141592653589793;

const int SEGCTION_LENGTH     = 10240;

const int NS_EQUATION         = 1;
const int TURBULENCE          = 2;

const int VIS_STD             = 0;
const int VIS_AVER            = 1;
const int VIS_TEST            = 2;
const int VIS_NEW1            = 3;
const int VIS_NEW2            = 4;
const int VIS_ERIC            = 5;

const int ISCHEME_ROE            = 1;
const int ISCHEME_VANLEER        = 2;
const int ISCHEME_STEGER         = 3;
const int ISCHEME_HLLE           = 4;
const int ISCHEME_LAX_FRIEDRICHS = 5;
const int ISCHEME_AUSM           = 6;
const int ISCHEME_AUSMP          = 7;
const int ISCHEME_AUSMPUP        = 8;
const int ISCHEME_AUSMDV         = 9;
const int ISCHEME_AUSM_W         = 10;
const int ISCHEME_AUSMPW         = 11;
const int ISCHEME_Rotate         = 12;
const int ISCHEME_ROE_MODIFIED   = 13;

const int FIRST_ORDER            = 1;
const int SECOND_ORDER           = 2;

const int ILMT_NOLIM             = 1;
const int ILMT_FIRST             = 2;
const int ILMT_BARTH             = 3;
const int ILMT_VENCAT            = 4;
const int ILMT_STRUCT            = 10;
const int ILMT_MINMOD            = 11;
const int ILMT_VAN_ALBADA        = 12;
const int ILMT_VAN_ALBADA_CLZ    = 13;
const int ILMT_VANLEER           = 14;
const int ILMT_SMOOTH            = 15;
const int ILMT_MINMODTEST        = 16;
const int ILMT_MIN_VAN           = 17;
const int ILMT_3rd_SMOOTH        = 18;
const int ILMT_3rd_Minmod_SMOOTH = 19;

const int INVISCID               = 0;
const int LAMINAR                = 1;
const int ALGEBRAIC              = 2;
const int ONE_EQU                = 3;
const int TWO_EQU                = 4;

const int BALDWIN_LOMAX          = 10;
const int SPALART_ALLMARAS       = 11;
const int KW_SST                 = 12;

const int MULTI_STAGE            = 1;
const int LU_SGS                 = 4;
const int Block_LU_SGS           = 5;

const int OneRefineOneDeform     = 1;
const int DeformAfterAllRefine   = 2;

const int X_DIR    = 1;
const int Y_DIR    = 2;
const int Z_DIR    = 3;

const int PHINT    = 1;
const int PHFLOAT  = 2;
const int PHDOUBLE = 3;
const int PHSTRING = 4;
const int PHBOOL   = 5;

const RDouble SAME_POINT_TOL = static_cast<RDouble> (1.0e-12);

namespace IDX
{
const int IR   = 0;
const int IU   = 1;
const int IV   = 2;
const int IW   = 3;
const int IP   = 4;

const int IRU  = 1;
const int IRV  = 2;
const int IRW  = 3;
const int IRE  = 4;

const int IA   = 0;
const int IX   = 1;
const int IY   = 2;
const int IZ   = 3;

const int IKES = 1;
const int IETA = 2;
const int IZTA = 3;

const int ISA  = 0;

const int IKE  = 0;
const int IKW  = 1;

const int ITT  = 0;
const int ITV  = 1;
const int ITE  = 2;

}
const int MAX_SPECIES_NUM = 20; //The maximum number of species.
const int MAX_REACTION_NUM = 50;//The maximum number of reactions.
}

namespace TEC_SPACE
{
const int ORDERED         = 0;
const int FELINESEG       = 1;
const int FETRIANGLE      = 2;
const int FEQUADRILATERAL = 3;
const int FETETRAHEDRON   = 4;
const int FEBRICK         = 5;
const int FEPOLYGON       = 6;
const int FEPOLYHEDRON    = 7;

const int WriteBoundary   = 0;
const int WriteBlock      = 1;
const int Writecomplete   = 2;

const int BoundaryVisual  = 0;
const int FieldVisual     = 1;
const int BlockVisual     = 2;
}