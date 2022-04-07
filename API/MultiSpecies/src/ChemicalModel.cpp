#include "ChemicalModel.h"
#include "Pointer.h"

namespace PHSPACE
{
    namespace GAS_SPACE
    {
        ChemicalModel::ChemicalModel()
        {
            nSpecies = 0;
            nReaction = 0;
            stoichForward = NULL;
            stoichBackward = NULL;
            collisionCoef = NULL;
            nReactionType = 0; //nType = 0 indicates the Gupta Model, nType = 1 indicates the Park Model, and nType = 2 indicates the Dunn-Kang Model.
            reactionData = NULL;
            reactionType = NULL;
        }

        ChemicalModel::ChemicalModel(int ns, int nr)
        {
            if (ns <= 0 || nr <= 0)
            {
                nSpecies = 0;
                nReaction = 0;
                stoichForward = NULL;
                stoichBackward = NULL;
                collisionCoef = NULL;
                nReactionType = 0;
                reactionData = NULL;
                reactionType = NULL;
            }
            else
            {
                nSpecies = ns;
                nReaction = nr;
                nReactionType = 0;
                AllocateMemories();
            }
        }

        ChemicalModel::~ChemicalModel()
        {
            DeAllocateMemories();
        }

        void ChemicalModel::AllocateMemories()
        {
            if (this->nSpecies <= 0 || this->nReaction <= 0)
                return;
            stoichForward = NewPointer2 <int> (this->nReaction, this->nSpecies);
            stoichBackward = NewPointer2 <int> (this->nReaction, this->nSpecies);
            collisionCoef = NewPointer2 <RDouble> (this->nReaction, this->nSpecies);
            reactionData = NewPointer2 <RDouble> (this->nReaction, 6);
            reactionType = NewPointer2 <int> (this->nReaction, 2);
        }

        void ChemicalModel::DeAllocateMemories()
        {
            if (this->nSpecies <= 0 || this->nReaction <= 0)
                return;
            // Remove the allocated memories.
            DelPointer2(stoichForward);
            DelPointer2(stoichBackward);
            DelPointer2(collisionCoef);
            DelPointer2(reactionData);
            DelPointer2(reactionType);

            nSpecies = 0;
            nReaction = 0;
            stoichForward = NULL;
            stoichBackward = NULL;
            collisionCoef = NULL;
            reactionData = NULL;
            reactionType = NULL;
        }

        void ChemicalModel::SetChemicalModelInfo(int ns, int nr)
        {
            if (ns <= 0 || nr <= 0)
                return;
            if (this->nSpecies == ns && this->nReaction == nr)
                return;
            DeAllocateMemories();

            this->nSpecies = ns;
            this->nReaction = nr;
            AllocateMemories();
        }

        void ChemicalModel::CreateGuptaModel(int ns/* = 11*/, int nr/* = 20*/)
        {
            this->SetChemicalReactionType(0);
            this->SetChemicalModelInfo(ns, nr);

            int *reactionCoef;
            RDouble *thirdBodyCoef, *rateCoef;
            if (ns == 5)
            {
                reactionCoef = new int[2 * 5 * 6];
                thirdBodyCoef = new RDouble [5 * 6];
                rateCoef = new RDouble [6 * 6];
                GenerateGupta5S6RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else if (ns == 7)
            {
                reactionCoef = new int[2 * 7 * 9];
                thirdBodyCoef = new RDouble [7 * 9];
                rateCoef = new RDouble [9 * 6];
                GenerateGupta7S9RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else if (ns == 11)
            {
                reactionCoef = new int[2 * 11 * 20];
                thirdBodyCoef = new RDouble [11 * 20];
                rateCoef = new RDouble [20 * 6];
                GenerateGupta11S20RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else
            {
                return;
            }

            int nIndex = 0, nStart = nr * ns;
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * ns;
                for (int s = 0; s < ns; ++ s)
                {
                    this->stoichForward[r][s]  = reactionCoef[nIndex + s];
                    this->stoichBackward[r][s] = reactionCoef[nStart + nIndex + s];
                    this->collisionCoef[r][s]  = thirdBodyCoef[nIndex + s];
                }
            }
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * 6;
                for (int s = 0; s < 6; ++ s)
                {
                    this->reactionData[r][s] = rateCoef[nIndex + s];
                }
            }

            //Set the type of collision reactions.
            for (int r = 0; r < nr; ++ r)
            {
                this->reactionType[r][0] = 0;
                for (int s = 0; s < ns; ++ s)
                {
                    if (this->collisionCoef[r][s] > 0.0)
                    {
                        this->reactionType[r][0] = 1;
                        break;
                    }
                }
            }

            delete []reactionCoef;
            delete []thirdBodyCoef;
            delete []rateCoef;
        }

        void ChemicalModel::CreateParkModel(int ns/* = 11*/, int nr/* = 20*/)
        {
            this->SetChemicalReactionType(1);
            this->SetChemicalModelInfo(ns, nr);

            int *reactionCoef;
            RDouble *thirdBodyCoef, *rateCoef;
            if (ns == 5)
            {
                reactionCoef = new int[2 * 5 * 17];
                thirdBodyCoef = new RDouble [5 * 17];
                rateCoef = new RDouble [17 * 6];
                GeneratePark5S17RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else if (ns == 7)
            {
                reactionCoef = new int[2 * 7 * 22];
                thirdBodyCoef = new RDouble [7 * 22];
                rateCoef = new RDouble [22 * 6];
                GeneratePark7S22RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else if (ns == 11)
            {
                reactionCoef = new int[2 * 11 * 48];
                thirdBodyCoef = new RDouble [11 * 48];
                rateCoef = new RDouble [48 * 6];
                GeneratePark11S48RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else
            {
                return;
            }

            int nIndex = 0, nStart = nr * ns;
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * ns;
                for (int s = 0; s < ns; ++ s)
                {
                    this->stoichForward[r][s]  = reactionCoef[nIndex + s];
                    this->stoichBackward[r][s] = reactionCoef[nStart + nIndex + s];
                    this->collisionCoef[r][s]  = thirdBodyCoef[nIndex + s];
                }
            }
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * 6;
                for (int s = 0; s < 6; ++ s)
                {
                    this->reactionData[r][s] = rateCoef[nIndex + s];
                }
            }

            //Set the type of collision reactions.
            for (int r = 0; r < nr; ++ r)
            {
                this->reactionType[r][0] = 0;
                for (int s = 0; s < ns; ++ s)
                {
                    if (this->collisionCoef[r][s] > 0.0)
                    {
                        this->reactionType[r][0] = 1;
                        break;
                    }
                }
            }

            delete []reactionCoef;
            delete []thirdBodyCoef;
            delete []rateCoef;
        }

        void ChemicalModel::CreateDunnKangModel(int ns/* = 11*/, int nr/* = 26*/)
        {
            this->SetChemicalReactionType(2);
            this->SetChemicalModelInfo(ns, nr);

            int *reactionCoef;
            RDouble *thirdBodyCoef, *rateCoef;
            if (ns == 5)
            {
                reactionCoef = new int[2 * 5 * 11];
                thirdBodyCoef = new RDouble [5 * 11];
                rateCoef = new RDouble [11 * 6];
                GenerateDunnKang5S11RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else if (ns == 7)
            {
                reactionCoef = new int[2 * 7 * 15];
                thirdBodyCoef = new RDouble [7 * 15];
                rateCoef = new RDouble [15 * 6];
                GenerateDunnKang7S15RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else if (ns == 11)
            {
                reactionCoef = new int[2 * 11 * 26];
                thirdBodyCoef = new RDouble [11 * 26];
                rateCoef = new RDouble [26 * 6];
                GenerateDunnKang11S26RData(reactionCoef, thirdBodyCoef, rateCoef);
            }
            else
            {
                return;
            }

            int nIndex = 0, nStart = nr * ns;
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * ns;
                for (int s = 0; s < ns; ++ s)
                {
                    this->stoichForward[r][s]  = reactionCoef[nIndex + s];
                    this->stoichBackward[r][s] = reactionCoef[nStart + nIndex + s];
                    this->collisionCoef[r][s]  = thirdBodyCoef[nIndex + s];
                }
            }
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * 6;
                for (int s = 0; s < 6; ++ s)
                {
                    this->reactionData[r][s] = rateCoef[nIndex + s];
                }
            }

            //Set the type of collision reactions.
            for (int r = 0; r < nr; ++ r)
            {
                this->reactionType[r][0] = 0;
                for (int s = 0; s < ns; ++ s)
                {
                    if (this->collisionCoef[r][s] > 0.0)
                    {
                        this->reactionType[r][0] = 1;
                        break;
                    }
                }
            }

            delete []reactionCoef;
            delete []thirdBodyCoef;
            delete []rateCoef;
        }

        void ChemicalModel::CreateMarsGasModel(int ns/* = 8*/, int nr/* = 12*/, int nType/* = 0*/)
        {
            this->SetChemicalReactionType(nType + 3);
            this->SetChemicalModelInfo(ns, nr);

            int *reactionCoef;
            RDouble *thirdBodyCoef, *rateCoef;
            if (ns == 8)
            {
                reactionCoef = new int[2 * 8 * 12];
                thirdBodyCoef = new RDouble [8 * 12];
                rateCoef = new RDouble [12 * 6];
                if (nType == 0)
                {
                    GenerateMarsPark8S12RData(reactionCoef, thirdBodyCoef, rateCoef);
                }
                else
                {
                    GenerateMarsMcKenzie8S12RData(reactionCoef, thirdBodyCoef, rateCoef);
                }
            }
            else
            {
                return;
            }

            int nIndex = 0, nStart = nr * ns;
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * ns;
                for (int s = 0; s < ns; ++ s)
                {
                    this->stoichForward[r][s]  = reactionCoef[nIndex + s];
                    this->stoichBackward[r][s] = reactionCoef[nStart + nIndex + s];
                    this->collisionCoef[r][s]  = thirdBodyCoef[nIndex + s];
                }
            }
            for (int r = 0; r < nr; ++ r)
            {
                nIndex = r * 6;
                for (int s = 0; s < 6; ++ s)
                {
                    this->reactionData[r][s] = rateCoef[nIndex + s];
                }
            }

            //Set the type of collision reactions.
            for (int r = 0; r < nr; ++ r)
            {
                this->reactionType[r][0] = 0;
                for (int s = 0; s < ns; ++ s)
                {
                    if (this->collisionCoef[r][s] > 0.0)
                    {
                        this->reactionType[r][0] = 1;
                        break;
                    }
                }
            }

            delete []reactionCoef;
            delete []thirdBodyCoef;
            delete []rateCoef;
        }

        void ChemicalModel::GenerateGupta5S6RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[6] = {0, 0, 0, 0, 2, 2};
            for (int r = 0; r < 6; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //!O, O2, NO, N, N2
            int stoichData[60] = {
                0,   1,   0,   0,   0,
                0,   0,   0,   0,   1,
                0,   0,   0,   1,   1,
                0,   0,   1,   0,   0,
                1,   0,   1,   0,   0,
                1,   0,   0,   0,   1,
                2,   0,   0,   0,   0,
                0,   0,   0,   2,   0,
                0,   0,   0,   3,   0,
                1,   0,   0,   1,   0,
                0,   1,   0,   1,   0,
                0,   0,   1,   1,   0 
            };
            RDouble collisionData[30] = {
                   25.0,   9.0,   1.0,   1.0,   2.0,
                    1.0,   1.0,   1.0,   0.0,   2.5,
                    0.0,   0.0,   0.0,   0.0,   0.0,
                   20.0,   1.0,  20.0,  20.0,   1.0,
                    0.0,   0.0,   0.0,   0.0,   0.0,
                    0.0,   0.0,   0.0,   0.0,   0.0 
            };
            RDouble rateData[36] = {
                3.6100E+12, -1.00,  5.9400E+04,  3.0100E+03, -0.50,  0.0000E+00,
                1.9200E+11, -0.50,  1.1310E+05,  1.0900E+04, -0.50,  0.0000E+00,
                4.1500E+16, -1.50,  1.1310E+05,  2.3200E+09, -1.50,  0.0000E+00,
                3.9700E+14, -1.50,  7.5600E+04,  1.0100E+08, -1.50,  0.0000E+00,
                3.1800E+03,  1.00,  1.9700E+04,  9.6300E+05,  0.50,  3.6000E+03,
                6.7500E+07,  0.00,  3.7500E+04,  1.5000E+07,  0.00,  0.0000E+00
            };

            for (int i = 0; i < 60; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 30; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 36; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GenerateGupta7S9RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[9] = {0, 0, 0, 0, 2, 2, 4, 4, 4};
            for (int r = 0; r < 9; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, NO+, N2, e-.
            int stoichData[126] = {
                0,   1,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   1,   0,
                0,   0,   1,   0,   0,   0,   0,
                1,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   1,   0,
                1,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   0,
                0,   0,   0,   3,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,
                0,   0,   1,   1,   0,   0,   0,
                0,   0,   0,   0,   1,   0,   1,
                0,   0,   1,   0,   1,   0,   1,
                0,   0,   0,   0,   1,   0,   1
            };
            RDouble collisionData[63] = {
                25.0,  9.0,   1.0,   1.0,   0.0,   2.0,   0.0,
                1.0,   1.0,   1.0,   0.0,   0.0,   2.5,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                20.0,  1.0,   20.0,  20.0,  0.0,   1.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   4.0,   0.0,   0.0,   0.0,   1.0,   0.0
            };
            RDouble rateData[54] = {
                3.6100E+12, -1.00,  5.9400E+04,  3.0100E+03, -0.50,  0.0000E+00,
                1.9200E+11, -0.50,  1.1310E+05,  1.0900E+04, -0.50,  0.0000E+00,
                4.1500E+16, -1.50,  1.1310E+05,  2.3200E+09, -1.50,  0.0000E+00,
                3.9700E+14, -1.50,  7.5600E+04,  1.0100E+08, -1.50,  0.0000E+00,
                3.1800E+03,  1.00,  1.9700E+04,  9.6300E+05,  0.50,  3.6000E+03,
                6.7500E+07,  0.00,  3.7500E+04,  1.5000E+07,  0.00,  0.0000E+00,
                9.0300E+03,  0.50,  3.2400E+04,  1.8000E+13, -1.00,  0.0000E+00,
                1.3800E+14, -1.84,  1.4100E+05,  1.0000E+12, -2.50,  0.0000E+00,
                2.2000E+09, -0.35,  1.0800E+05,  2.2000E+14, -2.50,  0.0000E+00
            };

            for (int i = 0; i < 126; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 63; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 54; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GenerateGupta11S20RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[20] = {0, 0, 0, 0, 2, 2, 4, 5, 5, 4, 3, 3, 4, 4, 4, 3, 3, 3, 3, 3};
            for (int r = 0; r < 20; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            int stoichData[440] = {
                0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,
                0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   1,
                2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   2,
                0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   2,
                0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,
                0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,
                0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   1,
                0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,
                0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,
                0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0
            };
            RDouble collisionData[220] = {
                25.0,  9.0,   1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   2.0,   0.0,
                1.0,   1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   2.5,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                20.0,  1.0,   20.0,  20.0,  0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   4.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0
            };
            RDouble rateData[120] = {
                3.6100E+12, -1.00,  5.9400E+04,  3.0100E+03, -0.50,  0.0000E+00,
                1.9200E+11, -0.50,  1.1310E+05,  1.0900E+04, -0.50,  0.0000E+00,
                4.1500E+16, -1.50,  1.1310E+05,  2.3200E+09, -1.50,  0.0000E+00,
                3.9700E+14, -1.50,  7.5600E+04,  1.0100E+08, -1.50,  0.0000E+00,
                3.1800E+03,  1.00,  1.9700E+04,  9.6300E+05,  0.50,  3.6000E+03,
                6.7500E+07,  0.00,  3.7500E+04,  1.5000E+07,  0.00,  0.0000E+00,
                9.0300E+03,  0.50,  3.2400E+04,  1.8000E+13, -1.00,  0.0000E+00,
                3.6000E+25, -2.91,  1.5800E+05,  2.2000E+28, -4.50,  0.0000E+00,
                1.1000E+26, -3.14,  1.6900E+05,  2.2000E+28, -4.50,  0.0000E+00,
                1.6000E+11, -0.98,  8.0800E+04,  8.0200E+15, -1.50,  0.0000E+00,
                2.9200E+12, -1.11,  2.8000E+04,  7.8000E+05,  0.50,  0.0000E+00,
                2.0200E+05,  0.81,  1.3000E+04,  7.8000E+05,  0.50,  0.0000E+00,
                1.4000E+07,  0.00,  6.7800E+04,  1.5000E+16, -1.50,  0.0000E+00,
                1.3800E+14, -1.84,  1.4100E+05,  1.0000E+12, -2.50,  0.0000E+00,
                2.2000E+09, -0.35,  1.0800E+05,  2.2000E+14, -2.50,  0.0000E+00,
                3.6300E+09, -0.60,  5.0800E+04,  1.5000E+07,  0.00,  0.0000E+00,
                3.4000E+13, -2.00,  2.3000E+04,  2.4800E+13, -2.20,  0.0000E+00,
                1.0000E+13, -0.93,  6.1000E+04,  4.8000E+08,  0.00,  0.0000E+00,
                1.8000E+09,  0.17,  3.3000E+04,  1.8000E+07,  0.50,  0.0000E+00,
                1.3400E+07,  0.31,  7.7270E+04,  1.0000E+08,  0.00,  0.0000E+00
            };

            for (int i = 0; i < 440; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 220; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 120; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GeneratePark5S17RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[17] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2};
            for (int r = 0; r < 17; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //!O, O2, NO, N, N2
            int stoichData[170] = {
                0,   2,   0,   0,   0,
                0,   1,   1,   0,   0,
                0,   1,   0,   0,   1,
                1,   1,   0,   0,   0,
                0,   1,   0,   1,   0,
                0,   1,   1,   0,   0,
                0,   0,   2,   0,   0,
                0,   0,   1,   0,   1,
                1,   0,   1,   0,   0,
                0,   0,   1,   1,   0,
                0,   1,   0,   0,   1,
                0,   0,   1,   0,   1,
                0,   0,   0,   0,   2,
                1,   0,   0,   0,   1,
                0,   0,   0,   1,   1,
                1,   0,   0,   0,   1,
                1,   0,   1,   0,   0,
                2,   1,   0,   0,   0,
                2,   0,   1,   0,   0,
                2,   0,   0,   0,   1,
                3,   0,   0,   0,   0,
                2,   0,   0,   1,   0,
                1,   1,   0,   1,   0,
                1,   0,   1,   1,   0,
                1,   0,   0,   1,   1,
                2,   0,   0,   1,   0,
                1,   0,   0,   2,   0,
                0,   1,   0,   2,   0,
                0,   0,   1,   2,   0,
                0,   0,   0,   2,   1,
                1,   0,   0,   2,   0,
                0,   0,   0,   3,   0,
                0,   0,   1,   1,   0,
                0,   1,   0,   1,   0
            };
            RDouble rateData[102] = {
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                6.4000E+11, -1.00,  3.8400E+04,  1.3159E+11, -1.00,  4.9744E+02,
                8.4000E+06,  0.00,  1.9450E+04,  3.2469E+07,  0.00,  3.4780E+03
            };

            for (int i = 0; i < 170; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 85; ++ i)
            {
                collisionCoef[i] = 0.0;
            }
            for (int i = 0; i < 102; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GeneratePark7S22RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[22] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 4, 1};
            for (int r = 0; r < 22; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, NO+, N2, e-.
            int stoichData[308] = {
                0,   2,   0,   0,   0,   0,   0,
                0,   1,   1,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   1,   0,
                1,   1,   0,   0,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   1,   0,   0,
                0,   1,   1,   0,   0,   0,   0,
                0,   0,   2,   0,   0,   0,   0,
                0,   0,   1,   0,   0,   1,   0,
                1,   0,   1,   0,   0,   0,   0,
                0,   0,   1,   1,   0,   0,   0,
                0,   0,   1,   0,   1,   0,   0,
                0,   1,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   2,   0,
                1,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   1,   0,
                0,   0,   0,   0,   1,   1,   0,
                1,   0,   0,   0,   0,   1,   0,
                1,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,
                0,   0,   0,   0,   0,   1,   1,
                2,   1,   0,   0,   0,   0,   0,
                2,   0,   1,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   1,   0,
                3,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   1,   0,   0,   0,
                2,   0,   0,   0,   1,   0,   0,
                1,   1,   0,   1,   0,   0,   0,
                1,   0,   1,   1,   0,   0,   0,
                1,   0,   0,   1,   0,   1,   0,
                2,   0,   0,   1,   0,   0,   0,
                1,   0,   0,   2,   0,   0,   0,
                1,   0,   0,   1,   1,   0,   0,
                0,   1,   0,   2,   0,   0,   0,
                0,   0,   1,   2,   0,   0,   0,
                0,   0,   0,   2,   0,   1,   0,
                1,   0,   0,   2,   0,   0,   0,
                0,   0,   0,   3,   0,   0,   0,
                0,   0,   0,   2,   1,   0,   0,
                0,   0,   1,   1,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,
                0,   0,   0,   0,   1,   0,   1,
                0,   0,   0,   2,   0,   0,   1
            };
            RDouble rateData[132] = {
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                6.4000E+11, -1.00,  3.8400E+04,  1.3159E+11, -1.00,  4.9744E+02,
                8.4000E+06,  0.00,  1.9450E+04,  3.2469E+07,  0.00,  3.4780E+03,
                8.8000E+02,  1.00,  3.3633E+04,  3.6997E+06,  1.00,  0.0000E+00,
                1.2000E+19, -1.60,  1.1392E+05,  1.2735E+08, -0.60,  0.0000E+00
            };

            for (int i = 0; i < 308; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 154; ++ i)
            {
                collisionCoef[i] = 0.0;
            }
            for (int i = 0; i < 132; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GeneratePark11S48RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[48] = {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,
                4, 4, 4, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1
            };
            for (int r = 0; r < 48; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            int stoichData[1056] = {
                0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   1,   0,   0,
                0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   0,   1,   0,
                1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,
                0,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,
                0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,
                0,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,
                0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   1,
                1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,
                0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,
                1,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,
                0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,
                2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                2,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,
                2,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,
                2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                1,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,
                2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,
                1,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,
                1,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0,
                0,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   2,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   1,   0,
                1,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   0,   1,   0,   0,   0,
                0,   0,   0,   2,   1,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,
                0,   0,   0,   2,   0,   1,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   1,   0,   0,   0,   0,
                0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,
                0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,
                0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   2,
                0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   2,
                0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   1,   0,   0,
                0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,
                0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   1
            };
            RDouble rateData[288] = {
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                1.0000E+16, -1.50,  6.0050E+04,  1.3353E+05, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                2.0000E+15, -1.50,  6.0050E+04,  2.6705E+04, -0.50,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                1.1000E+11,  0.00,  7.6022E+04,  5.6774E+00,  1.00,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                5.0000E+09,  0.00,  7.6022E+04,  2.5807E-01,  1.00,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                3.0000E+16, -1.60,  1.1392E+05,  3.1837E+05, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                7.0000E+15, -1.60,  1.1392E+05,  7.4287E+04, -0.60,  0.0000E+00,
                6.4000E+11, -1.00,  3.8400E+04,  1.3159E+11, -1.00,  4.9744E+02,
                8.4000E+06,  0.00,  1.9450E+04,  3.2469E+07,  0.00,  3.4780E+03,
                8.8000E+02,  1.00,  3.3633E+04,  3.6997E+06,  1.00,  0.0000E+00,
                4.4000E+01,  1.50,  6.9355E+04,  4.8044E+04,  1.50,  0.0000E+00,
                7.1000E-04,  2.70,  8.2471E+04,  2.1358E+00,  2.70,  0.0000E+00,
                3.9000E+27, -3.78,  1.6108E+05,  2.0798E+20, -2.78,  0.0000E+00,
                2.5000E+28, -3.82,  1.7175E+05,  1.0557E+21, -2.82,  0.0000E+00,
                1.0000E+06,  0.50,  7.8067E+04,  2.3385E+05,  0.50,  0.0000E+00,
                8.7000E+07,  0.14,  2.9229E+04,  2.8434E+07,  0.14,  0.0000E+00,
                1.4000E-01,  1.90,  2.6644E+04,  1.3322E-01,  1.90,  0.0000E+00,
                9.9000E+06,  0.00,  4.0759E+04,  2.8561E+06,  0.00,  0.0000E+00,
                4.0000E+06, -0.09,  1.8557E+04,  5.3106E+06, -0.09,  0.0000E+00,
                3.4000E+07, -1.08,  1.3520E+04,  4.0638E+07, -1.08,  0.0000E+00,
                2.4000E+07,  0.41,  3.2866E+04,  4.4426E+06,  0.41,  0.0000E+00,
                7.2000E+06,  0.29,  4.8838E+04,  5.1518E+06,  0.29,  0.0000E+00,
                9.1000E+05,  0.36,  2.2800E+04,  1.9774E+05,  0.36,  5.9851E+02,
                7.2000E+07,  0.00,  3.5722E+04,  1.8700E+07,  0.00,  0.0000E+00,
                1.2000E+19, -1.60,  1.1392E+05,  1.2735E+08, -0.60,  0.0000E+00
            };

            for (int i = 0; i < 1056; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 528; ++ i)
            {
                collisionCoef[i] = 0.0;
            }
            for (int i = 0; i < 288; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GenerateDunnKang5S11RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2};
            for (int r = 0; r < 11; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //!O, O2, NO, N, N2
            int stoichData[110] = {
                0,   1,   0,   0,   0,
                1,   1,   0,   0,   0,
                0,   2,   0,   0,   0,
                0,   1,   0,   0,   1,
                0,   0,   0,   0,   1,
                0,   0,   0,   1,   1,
                0,   0,   0,   0,   2,
                0,   0,   1,   0,   0,
                0,   0,   1,   0,   0,
                1,   0,   1,   0,   0,
                1,   0,   0,   0,   1,
                2,   0,   0,   0,   0,
                3,   0,   0,   0,   0,
                2,   1,   0,   0,   0,
                2,   0,   0,   0,   1,
                0,   0,   0,   2,   0,
                0,   0,   0,   3,   0,
                0,   0,   0,   2,   1,
                1,   0,   0,   1,   0,
                1,   0,   0,   1,   0,
                0,   1,   0,   1,   0,
                0,   0,   1,   1,   0
            };
            RDouble collisionData[55] = {
                0,   0,   1,   1,   0,
                0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,
                1,   1,   1,   0,   0,
                0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,
                0,   1,   0,   0,   1,
                1,   0,   1,   1,   0,
                0,   0,   0,   0,   0,
                0,   0,   0,   0,   0
            };
            RDouble rateData[66] = {
                3.6000E+12, -1.00,  5.9500E+04,  3.0000E+03, -0.50,  0.0000E+00,
                9.0000E+13, -1.00,  5.9500E+04,  7.5000E+04, -0.50,  0.0000E+00,
                3.2400E+13, -1.00,  5.9500E+04,  2.7000E+04, -0.50,  0.0000E+00,
                7.2000E+12, -1.00,  5.9500E+04,  6.0000E+03, -0.50,  0.0000E+00,
                1.9000E+11, -0.50,  1.1300E+05,  1.1000E+04, -0.50,  0.0000E+00,
                4.0850E+16, -1.50,  1.1300E+05,  2.2700E+09, -1.50,  0.0000E+00,
                4.7000E+11, -0.50,  1.1300E+05,  2.7200E+04, -0.50,  0.0000E+00,
                3.9000E+14, -1.50,  7.5500E+04,  1.0000E+08, -1.50,  0.0000E+00,
                7.8000E+14, -1.50,  7.5500E+04,  2.0000E+08, -1.50,  0.0000E+00,
                3.2000E+03,  1.00,  1.9700E+04,  1.3000E+04,  1.00,  3.5800E+03,
                7.0000E+07,  0.00,  3.8000E+04,  1.5600E+07,  0.00,  0.0000E+00
            };

            for (int i = 0; i < 110; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 55; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 66; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GenerateDunnKang7S15RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 4, 4, 4, 4};
            for (int r = 0; r < 15; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, NO+, N2, e-.
            int stoichData[210] = {
                0,   1,   0,   0,   0,   0,   0,
                1,   1,   0,   0,   0,   0,   0,
                0,   2,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   1,   0,
                0,   0,   0,   0,   0,   2,   0,
                0,   0,   1,   0,   0,   0,   0,
                0,   0,   1,   0,   0,   0,   0,
                1,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   1,   0,
                1,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   1,   0,
                0,   1,   1,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,
                3,   0,   0,   0,   0,   0,   0,
                2,   1,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   2,   0,   0,   0,
                0,   0,   0,   3,   0,   0,   0,
                0,   0,   0,   2,   0,   1,   0,
                1,   0,   0,   1,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,
                0,   0,   1,   1,   0,   0,   0,
                0,   0,   0,   0,   1,   0,   1,
                0,   0,   1,   0,   1,   0,   1,
                0,   0,   0,   0,   1,   1,   1,
                0,   1,   0,   0,   1,   0,   1
            };
            RDouble collisionData[105] = {
                0,   0,   1,   1,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                1,   1,   1,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   1,   0,
                1,   0,   1,   1,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0
            };
            RDouble rateData[90] = {
                3.6000E+12, -1.00,  5.9500E+04,  3.0000E+03, -0.50,  0.0000E+00,
                9.0000E+13, -1.00,  5.9500E+04,  7.5000E+04, -0.50,  0.0000E+00,
                3.2400E+13, -1.00,  5.9500E+04,  2.7000E+04, -0.50,  0.0000E+00,
                7.2000E+12, -1.00,  5.9500E+04,  6.0000E+03, -0.50,  0.0000E+00,
                1.9000E+11, -0.50,  1.1300E+05,  1.1000E+04, -0.50,  0.0000E+00,
                4.0850E+16, -1.50,  1.1300E+05,  2.2700E+09, -1.50,  0.0000E+00,
                4.7000E+11, -0.50,  1.1300E+05,  2.7200E+04, -0.50,  0.0000E+00,
                3.9000E+14, -1.50,  7.5500E+04,  1.0000E+08, -1.50,  0.0000E+00,
                7.8000E+14, -1.50,  7.5500E+04,  2.0000E+08, -1.50,  0.0000E+00,
                3.2000E+03,  1.00,  1.9700E+04,  1.3000E+04,  1.00,  3.5800E+03,
                7.0000E+07,  0.00,  3.8000E+04,  1.5600E+07,  0.00,  0.0000E+00,
                1.4000E+00,  1.50,  3.1900E+04,  6.7000E+15, -1.50,  0.0000E+00,
                1.3800E+14, -1.84,  1.4100E+05,  1.0000E+12, -2.50,  0.0000E+00,
                2.2000E+09, -0.35,  1.0800E+05,  2.2000E+14, -2.50,  0.0000E+00,
                8.8000E+09, -0.35,  1.0800E+05,  8.8000E+14, -2.50,  0.0000E+00
            };

            for (int i = 0; i < 210; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 105; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 90; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GenerateDunnKang11S26RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[26] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 4, 4, 4, 4, 5, 5, 3, 3, 3, 3, 3, 4, 4, 3, 3};
            for (int r = 0; r < 26; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            int stoichData[572] = {
                0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   1,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,
                0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   1,
                0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                2,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   2,   0,   0,   0,   0,   0,   1,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,
                0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   1,
                0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   1,
                0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   1,
                0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   2,
                0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   2,
                0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,
                0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,
                1,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,
                0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,
                0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,
                0,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0
            };
            RDouble collisionData[286] = {
                0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
                1,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
            };
            RDouble rateData[156] = {
                3.6000E+12, -1.00,  5.9500E+04,  3.0000E+03, -0.50,  0.0000E+00,
                9.0000E+13, -1.00,  5.9500E+04,  7.5000E+04, -0.50,  0.0000E+00,
                3.2400E+13, -1.00,  5.9500E+04,  2.7000E+04, -0.50,  0.0000E+00,
                7.2000E+12, -1.00,  5.9500E+04,  6.0000E+03, -0.50,  0.0000E+00,
                1.9000E+11, -0.50,  1.1300E+05,  1.1000E+04, -0.50,  0.0000E+00,
                4.0850E+16, -1.50,  1.1300E+05,  2.2700E+09, -1.50,  0.0000E+00,
                4.7000E+11, -0.50,  1.1300E+05,  2.7200E+04, -0.50,  0.0000E+00,
                3.9000E+14, -1.50,  7.5500E+04,  1.0000E+08, -1.50,  0.0000E+00,
                7.8000E+14, -1.50,  7.5500E+04,  2.0000E+08, -1.50,  0.0000E+00,
                3.2000E+03,  1.00,  1.9700E+04,  1.3000E+04,  1.00,  3.5800E+03,
                7.0000E+07,  0.00,  3.8000E+04,  1.5600E+07,  0.00,  0.0000E+00,
                1.4000E+00,  1.50,  3.1900E+04,  6.7000E+15, -1.50,  0.0000E+00,
                1.3800E+14, -1.84,  1.4100E+05,  1.0000E+12, -2.50,  0.0000E+00,
                2.2000E+09, -0.35,  1.0800E+05,  2.2000E+14, -2.50,  0.0000E+00,
                8.8000E+09, -0.35,  1.0800E+05,  8.8000E+14, -2.50,  0.0000E+00,
                3.6000E+25, -2.91,  1.5800E+05,  2.2000E+28, -4.50,  0.0000E+00,
                1.1000E+26, -3.14,  1.6900E+05,  2.2000E+28, -4.50,  0.0000E+00,
                1.0000E+13, -0.93,  6.1000E+04,  4.8000E+08,  0.00,  0.0000E+00,
                1.3400E+07,  0.31,  7.7270E+04,  1.0000E+08,  0.00,  0.0000E+00,
                3.6300E+09, -0.60,  5.0800E+04,  1.5000E+07,  0.00,  0.0000E+00,
                3.4000E+13, -2.00,  2.3000E+04,  2.4800E+13, -2.20,  0.0000E+00,
                2.0200E+05,  0.81,  1.3000E+04,  7.8000E+05,  0.50,  0.0000E+00,
                1.4000E+07,  0.00,  6.7800E+04,  1.5000E+16, -1.50,  0.0000E+00,
                1.6000E+11, -0.98,  8.0800E+04,  8.0000E+15, -1.50,  0.0000E+00,
                1.8000E+09,  0.17,  3.3000E+04,  1.8000E+07,  0.50,  0.0000E+00,
                2.9200E+12, -1.11,  2.8000E+04,  7.8000E+05,  0.50,  0.0000E+00
            };

            for (int i = 0; i < 572; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 286; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 156; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GenerateMarsPark8S12RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[12] = {0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2};
            for (int r = 0; r < 12; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, N2, C, CO, CO2.
            int stoichData[192] = {
                0,  0,  0,  0,  0,  0,  0,  1,
                0,  0,  0,  0,  0,  0,  1,  0,
                0,  0,  0,  0,  1,  0,  0,  0,
                0,  1,  0,  0,  0,  0,  0,  0,
                0,  0,  1,  0,  0,  0,  0,  0,
                0,  0,  0,  1,  0,  0,  1,  0,
                0,  0,  0,  0,  0,  0,  2,  0,
                1,  0,  0,  0,  0,  0,  1,  0,
                1,  0,  0,  0,  1,  0,  0,  0,
                0,  0,  1,  0,  0,  0,  1,  0,
                1,  0,  1,  0,  0,  0,  0,  0,
                1,  0,  0,  0,  0,  0,  0,  1,
                1,  0,  0,  0,  0,  0,  1,  0,
                1,  0,  0,  0,  0,  1,  0,  0,
                0,  0,  0,  2,  0,  0,  0,  0,
                2,  0,  0,  0,  0,  0,  0,  0,
                1,  0,  0,  1,  0,  0,  0,  0,
                0,  0,  1,  0,  0,  1,  0,  0,
                0,  0,  0,  0,  0,  1,  0,  1,
                0,  1,  0,  0,  0,  1,  0,  0,
                0,  0,  1,  1,  0,  0,  0,  0,
                0,  0,  0,  1,  0,  0,  0,  1,
                0,  1,  0,  1,  0,  0,  0,  0,
                0,  1,  0,  0,  0,  0,  1,  0
            };
            RDouble collisionData[96] = {
                1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,
                14.80, 1.00,  1.00,  14.80, 1.00,  14.80, 1.95,  1.00,
                4.28,  1.00,  1.00,  4.28,  1.00,  4.28,  1.00,  1.00,
                5.00,  1.00,  1.00,  5.00,  1.00,  5.00,  1.00,  1.00,
                22.00, 1.00,  22.00, 22.00, 1.00,  22.00, 1.00,  22.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00
            };
            RDouble rateData[72] = {
                3.70E+08,   0.00,  5.250E+04,  6.10E-03,   0.75,  -1.024E+04,
                2.30E+13,  -1.00,  1.290E+05,  5.13E+05,  -1.00,   0.000E+00,
                7.00E+15,  -1.60,  1.132E+05,  4.27E+08,  -1.60,   0.000E+00,
                2.00E+15,  -1.50,  5.950E+04,  2.00E+06,  -1.00,   0.000E+00,
                5.00E+09,   0.00,  7.550E+04,  4.27E+03,   0.00,   0.000E+00,
                2.86E+05,   0.50,  5.363E+04,  2.60E+04,   0.50,   0.000E+00,
                2.33E+03,   0.50,  6.571E+04,  4.60E+06,  -0.25,   0.000E+00,
                3.90E+07,  -0.18,  6.920E+04,  1.34E+08,  -0.43,   0.000E+00,
                6.40E+11,  -1.00,  3.830E+04,  1.39E+11,  -1.00,   0.000E+00,
                4.59E+02,   0.50,  1.207E+04,  9.90E+06,  -0.25,   0.000E+00,
                8.40E+06,   0.00,  1.945E+04,  2.68E+05,   0.50,   0.000E+00,
                1.70E+07,   0.00,  2.650E+04,  3.33E+04,   0.50,   2.250E+04
            };

            for (int i = 0; i < 192; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 96; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 72; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        void ChemicalModel::GenerateMarsMcKenzie8S12RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef)
        {
            //Set the types of chemical reactions.
            int typeOfReaction[12] = {0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2};
            for (int r = 0; r < 12; ++ r)
            {
                reactionType[r][1] = typeOfReaction[r];
            }

            //! O, O2, NO, N, N2, C, CO, CO2.
            int stoichData[192] = {
                0,  0,  0,  0,  0,  0,  0,  1,
                0,  0,  0,  0,  0,  0,  1,  0,
                0,  0,  0,  0,  1,  0,  0,  0,
                0,  1,  0,  0,  0,  0,  0,  0,
                0,  0,  1,  0,  0,  0,  0,  0,
                0,  0,  0,  1,  0,  0,  1,  0,
                0,  0,  0,  0,  0,  0,  2,  0,
                1,  0,  0,  0,  0,  0,  1,  0,
                1,  0,  0,  0,  1,  0,  0,  0,
                0,  0,  1,  0,  0,  0,  1,  0,
                1,  0,  1,  0,  0,  0,  0,  0,
                1,  0,  0,  0,  0,  0,  0,  1,
                1,  0,  0,  0,  0,  0,  1,  0,
                1,  0,  0,  0,  0,  1,  0,  0,
                0,  0,  0,  2,  0,  0,  0,  0,
                2,  0,  0,  0,  0,  0,  0,  0,
                1,  0,  0,  1,  0,  0,  0,  0,
                0,  0,  1,  0,  0,  1,  0,  0,
                0,  0,  0,  0,  0,  1,  0,  1,
                0,  1,  0,  0,  0,  1,  0,  0,
                0,  0,  1,  1,  0,  0,  0,  0,
                0,  0,  0,  1,  0,  0,  0,  1,
                0,  1,  0,  1,  0,  0,  0,  0,
                0,  1,  0,  0,  0,  0,  1,  0
            };
            RDouble collisionData[96] = {
                1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,
                14.80, 1.00,  1.00,  14.80, 1.00,  14.80, 1.95,  1.00,
                4.28,  1.00,  1.00,  4.28,  1.00,  4.28,  1.00,  1.00,
                5.00,  1.00,  1.00,  5.00,  1.00,  5.00,  1.00,  1.00,
                22.00, 1.00,  22.00, 22.00, 1.00,  22.00, 1.00,  22.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00
            };
            RDouble rateData[72] = {
                1.20E+05,   0.50,   3.434E+04,   2.00E-06,   1.25,  -2.840E+04,
                4.48E+13,  -1.00,   1.289E+05,   1.00E+06,  -1.00,   0.000E+00,
                2.46E+13,  -1.00,   1.132E+05,   1.50E+06,  -1.00,   0.000E+00,
                9.05E+12,  -1.00,   5.937E+04,   9.00E+03,  -0.50,   0.000E+00,
                4.09E+12,  -1.00,   7.533E+04,   3.50E+06,  -1.00,   0.000E+00,
                2.86E+05,   0.50,   5.360E+04,   2.60E+04,   0.50,   0.000E+00,
                2.33E+03,   0.50,   6.571E+04,   4.60E+06,  -0.25,   0.000E+00,
                2.73E+06,   0.50,   6.954E+04,   9.40E+06,   0.25,   0.000E+00,
                7.35E+05,   0.50,   3.794E+04,   1.60E+05,   0.50,   0.000E+00,
                4.59E+02,   0.50,   1.207E+04,   9.90E+06,  -0.25,   0.000E+00,
                2.98E+05,   0.50,   1.946E+04,   9.50E+03,   1.00,   0.000E+00,
                2.54E+05,   0.50,   2.769E+04,   5.00E+02,   1.00,   2.390E+04,
            };

            for (int i = 0; i < 192; ++ i)
            {
                reactionCoef[i] = stoichData[i];
            }
            for (int i = 0; i < 96; ++ i)
            {
                collisionCoef[i] = collisionData[i];
            }
            for (int i = 0; i < 72; ++ i)
            {
                rateCoef[i] = rateData[i];
            }
        }

        //===============================================
        ChemkinData::ChemkinData()
        {
            nSpecies = 0;//number of species;
            nInterval = 0;
            nCoefNum = 0;
            temperatureRange = NULL;
            curveFittingData = NULL;
        }

        ChemkinData::~ChemkinData()
        {
            DeAllocateMemories();
        }

        void ChemkinData::AllocateMemories()
        {
            if (nSpecies > 0 && nInterval > 0 && nCoefNum > 0)
            {
                temperatureRange = new RDouble[nInterval + 1];
                curveFittingData = NewPointer3 <RDouble> (this->nSpecies, nInterval, nCoefNum);
            }
        }

        void ChemkinData::DeAllocateMemories()
        {
            if (nSpecies > 0 && nInterval > 0 && nCoefNum > 0)
            {
                delete []temperatureRange;
                DelPointer3(curveFittingData);
                nSpecies = 0;//number of species;
                curveFittingData = NULL;
                nInterval = 0;
                nCoefNum = 0;
                temperatureRange = NULL;
            }
        }

        void ChemkinData::SetInterpolationInfo(int ns, int m/* = 5*/, int n/* = 6*/)
        {
            if (ns <= 0 || m <= 0 || n <= 0)
                return;
            if (ns == this->nSpecies && m == this->nInterval && n == this->nCoefNum)
                return;
            DeAllocateMemories();

            this->nSpecies = ns;
            this->nInterval = m;
            this->nCoefNum  = n;
            AllocateMemories();
        }

        void ChemkinData::CreateChemkinData(int nGasModel/* = 0*/, int ns/* = 11*/)
        {
            RDouble ***fittingData;
            this->SetInterpolationInfo(ns, 5, 6);
            RDouble temperature[6] = {300.0, 1000.0, 6000.0, 15000.0, 25000.0, 30000.0};
            for (int i = 0; i < 6; ++ i)
            {
                temperatureRange[i] = temperature[i];
            }
            fittingData = NewPointer3 <RDouble> (16, 5, 6);
            // Generate the original data.
            this->GenerateFullChemkinData(fittingData);

            if (nGasModel == 0 && (ns == 5 || ns == 7)) //The earth atmosphere.
            {
                // To create the curve fitting data of specific model.
                for (int s = 0; s < 4; ++ s)
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[s][i][j] = fittingData[s][i][j];
                        }
                    }
                }
                if (ns == 5)
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[4][i][j] = fittingData[9][i][j];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[4][i][j] = fittingData[6][i][j];
                            this->curveFittingData[5][i][j] = fittingData[9][i][j];
                            this->curveFittingData[6][i][j] = fittingData[10][i][j];
                        }
                    }
                }
            }
            else if (nGasModel == 0 && ns == 11) //The earth atmosphere.
            {
                for (int s = 0; s < 11; ++ s)
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[s][i][j] = fittingData[s][i][j];
                        }
                    }
                }
            }
            else if (nGasModel == 1 && ns == 8) //The mars atmosphere.
            {
                for (int s = 0; s < 4; ++ s)
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[s][i][j] = fittingData[s][i][j];
                        }
                    }
                }
                for (int i = 0; i < 5; ++ i)
                {
                    for (int j = 0; j < 6; ++ j)
                    {
                        this->curveFittingData[4][i][j] = fittingData[9][i][j];
                        this->curveFittingData[5][i][j] = fittingData[11][i][j];
                        this->curveFittingData[6][i][j] = fittingData[14][i][j];
                        this->curveFittingData[7][i][j] = fittingData[15][i][j];
                    }
                }
            }
            else if (nGasModel == 1 && ns == 10) //The mars atmosphere.
            {
                for (int s = 0; s < 4; ++ s)
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[s][i][j] = fittingData[s][i][j];
                        }
                    }
                }
                for (int i = 0; i < 5; ++ i)
                {
                    for (int j = 0; j < 6; ++ j)
                    {
                        this->curveFittingData[4][i][j] = fittingData[9][i][j];
                    }
                }
                for (int s = 11; s < 16; ++ s)
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[s - 6][i][j] = fittingData[s][i][j];
                        }
                    }
                }
            }
            else
            {
                for (int s = 0; s < ns; ++ s)
                {
                    for (int i = 0; i < 5; ++ i)
                    {
                        for (int j = 0; j < 6; ++ j)
                        {
                            this->curveFittingData[s][i][j] = fittingData[s][i][j];
                        }
                    }
                }
            }
            DelPointer3(fittingData);
        }

        void ChemkinData::GenerateFullChemkinData(RDouble ***fittingData)
        {
            // the original data.
            // the order of species is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-, C, C2, CN, CO, CO2.
            RDouble curveFitData[480] = {
                2.823600E+00, -8.947800E-04,  8.306000E-07, -1.683700E-10, -7.320500E-14,  2.915000E+04,
                2.542100E+00, -2.755100E-05, -3.102800E-09,  4.551100E-12, -4.368100E-16,  2.915000E+04,
                2.546000E+00, -5.952000E-05,  2.701000E-08, -2.798000E-12,  9.380000E-17,  2.915000E+04,
               -9.787100E-03,  1.245000E-03, -1.615400E-07,  8.038000E-12, -1.262400E-16,  2.915000E+04,
                1.642800E+01, -3.931300E-03,  2.984000E-07, -8.161300E-12,  7.500400E-17,  2.915000E+04,

                3.614600E+00, -1.859800E-03,  7.081400E-06, -6.807000E-09,  2.162800E-12, -1.044000E+03,
                3.594900E+00,  7.521300E-04, -1.873200E-07,  2.791300E-11, -1.577400E-15, -1.044000E+03,
                3.859900E+00,  3.251000E-04, -9.213100E-09, -7.868400E-13,  2.942600E-17, -1.044000E+03,
                3.486700E+00,  5.238400E-04, -3.912300E-08,  1.009400E-12, -8.871800E-18, -1.044000E+03,
                3.962000E+00,  3.944600E-04, -2.950600E-08,  7.397500E-13, -6.420900E-18, -1.044000E+03,

                3.588700E+00, -1.247900E-03,  3.978600E-06, -2.865100E-09,  6.301500E-13,  9.764000E+03,
                3.204700E+00,  1.270500E-03, -4.660300E-07,  7.500700E-11, -4.231400E-15,  9.764000E+03,
                3.854300E+00,  2.340900E-04, -2.135400E-08,  1.668900E-12, -4.907000E-17,  9.764000E+03,
                4.330900E+00, -5.808600E-05,  2.805900E-08, -1.569400E-12,  2.410400E-17,  9.764000E+03,
                2.350700E+00,  5.864300E-04, -3.131600E-08,  6.049500E-13, -4.055700E-18,  9.764000E+03,

                2.503100E+00, -2.180000E-05,  5.420500E-08, -5.647600E-11,  2.099900E-14,  5.613000E+04,
                2.482000E+00,  6.925800E-05, -6.306500E-08,  1.838700E-11, -1.174700E-15,  5.613000E+04,
                2.748000E+00, -3.909000E-04,  1.338000E-07, -1.191000E-11,  3.369000E-16,  5.613000E+04,
               -1.228000E+00,  1.926800E-03, -2.437000E-07,  1.219300E-11, -1.991800E-16,  5.613000E+04,
                1.552000E+01, -3.885800E-03,  3.228800E-07, -9.605300E-12,  9.547200E-17,  5.613000E+04,

                2.498500E+00,  1.141100E-05, -2.976100E-08,  3.224700E-11, -1.237600E-14,  1.879000E+05,
                2.506000E+00, -1.446400E-05,  1.244600E-08, -4.685800E-12,  6.554900E-16,  1.879000E+05,
                2.944000E+00, -4.108000E-04,  9.156000E-08, -5.848000E-12,  1.190000E-16,  1.879000E+05,
                1.278400E+00,  4.086600E-04, -2.173100E-08,  3.325200E-13,  6.316000E-19,  1.879000E+05,
                1.288900E+00,  4.334300E-04, -2.675800E-08,  6.215900E-13, -4.513100E-18,  1.879000E+05,

                3.243000E+00,  1.174000E-03, -3.900000E-07,  5.437000E-11, -2.392000E-15,  1.400000E+05,
                3.243000E+00,  1.174000E-03, -3.900000E-07,  5.437000E-11, -2.392000E-15,  1.400000E+05,
                5.169000E+00, -8.620000E-04,  2.041000E-07, -1.300000E-11,  2.494000E-16,  1.400000E+05,
               -2.801700E-01,  1.667400E-03, -1.210700E-07,  3.211300E-12, -2.834900E-17,  1.400000E+05,
                2.044500E+00,  1.031300E-03, -7.404600E-08,  1.925700E-12, -1.746100E-17,  1.400000E+05,

                3.529400E+00, -3.034200E-04,  3.854400E-07,  1.051900E-09, -7.277700E-13,  1.184000E+05,
                3.215200E+00,  9.974200E-04, -2.903000E-07,  3.692500E-11, -1.599400E-15,  1.184000E+05,
                2.689600E+00,  1.379600E-03, -3.398500E-07,  3.377600E-11, -1.042700E-15,  1.184000E+05,
                5.934600E+00, -1.317800E-03,  2.329700E-07, -1.173300E-11,  1.840200E-16,  1.184000E+05,
               -5.159500E+00,  2.629000E-03, -1.625400E-07,  3.938100E-12, -3.431100E-17,  1.184000E+05,

                2.727000E+00, -2.820000E-04,  1.105000E-07, -1.551000E-11,  7.847000E-16,  2.254000E+05,
                2.727000E+00, -2.820000E-04,  1.105000E-07, -1.551000E-11,  7.847000E-16,  2.254000E+05,
                2.499000E+00, -3.725000E-06,  1.147000E-08, -1.102000E-12,  3.078000E-17,  2.254000E+05,
                2.385600E+00,  8.349500E-05, -5.881500E-09,  1.885000E-13, -1.612000E-18,  2.254000E+05,
                2.228600E+00,  1.245800E-04, -8.763600E-09,  2.620400E-13, -2.167400E-18,  2.254000E+05,

                3.549800E+00, -6.081000E-04,  1.469000E-06, -6.509100E-11, -3.564900E-13,  1.826000E+05,
                3.397000E+00,  4.525000E-04,  1.272000E-07, -3.879000E-11,  2.459000E-15,  1.826000E+05,
                3.378000E+00,  8.629000E-04, -1.276000E-07,  8.087000E-12, -1.880000E-16,  1.826000E+05,
                4.394200E+00,  1.886800E-04, -7.127200E-09, -1.751100E-13,  6.717600E-18,  1.826000E+05,
                3.949300E+00,  3.679500E-04, -2.691000E-08,  6.711000E-13, -5.824400E-18,  1.826000E+05,

                3.674800E+00, -1.208100E-03,  2.324000E-06, -6.321800E-10, -2.257700E-13, -1.043000E+03,
                3.212500E+00,  1.013700E-03, -3.046700E-07,  4.109100E-11, -2.017000E-15, -1.043000E+03,
                3.181100E+00,  8.974500E-04, -2.021600E-07,  1.826600E-11, -5.033400E-16, -1.043000E+03,
                9.637700E+00, -2.572800E-03,  3.302000E-07, -1.431500E-11,  2.033300E-16, -1.043000E+03,
               -5.168100E+00,  2.333700E-03, -1.295300E-07,  2.787200E-12, -2.136000E-17, -1.043000E+03,

                2.500000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00, -7.454200E+02,
                2.500000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00, -7.454200E+02,
                2.500000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00, -7.454200E+02,
                2.500000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00, -7.454200E+02,
                2.500000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00, -7.454200E+02,

                2.50000000E+00,  2.91840000E-08,  0.00000000E+00,  0.00000000E+00,  0.00000000E+00,  8.62400000E+04,
                2.73010000E+00, -3.54680000E-04,  1.76380000E-07, -2.90810000E-11,  1.61410000E-15,  8.62400000E+04,
                2.19260000E+00,  2.33960000E-04, -3.29930000E-08,  1.86950000E-12, -3.85470000E-17,  8.62400000E+04,
                3.09520000E+00, -5.78480000E-05,  2.53490000E-09, -5.48800000E-14,  4.76900000E-19,  8.62400000E+04,
                2.99920000E+00, -4.21660000E-05,  1.57100000E-09, -2.84690000E-14,  2.04760000E-19,  8.62400000E+04,

                3.59270000E+00, -1.31500000E-03,  4.68360000E-06, -3.64550000E-09,  8.22490000E-13,  9.88690000E+04,
                3.46870000E+00,  9.91340000E-04, -3.80030000E-07,  6.51940000E-11, -3.95530000E-15,  9.88690000E+04,
                3.70860000E+00,  2.26880000E-04, -9.33260000E-09, -3.05890000E-13,  1.81940000E-17,  9.88690000E+04,
                4.40650000E+00,  1.37050000E-04, -1.14540000E-08,  3.62640000E-13, -4.12340000E-18,  9.88690000E+04,
                5.45660000E+00, -4.88060000E-05,  9.16120000E-10, -4.04850000E-15, -4.16240000E-20,  9.88690000E+04,

                3.52730000E+00, -2.94470000E-04,  3.54480000E-07,  1.61510000E-09, -1.20290000E-12,  5.25620000E+04,
                3.76220000E+00,  3.96600000E-05,  2.57790000E-07, -6.11840000E-11,  3.96340000E-15,  5.25620000E+04,
                5.08050000E+00,  2.03760000E-04, -4.76650000E-08,  3.31230000E-12, -7.70480000E-17,  5.25620000E+04,
                6.03680000E+00, -1.99970000E-04,  1.09990000E-08, -2.86890000E-13,  2.91880000E-18,  5.25620000E+04,
                5.36840000E+00, -8.39690000E-05,  3.41690000E-09, -6.58850000E-14,  4.96200000E-19,  5.25620000E+04,

                3.51580000E+00, -1.30220000E-04, -2.62670000E-07,  2.25550000E-09, -1.41070000E-12, -1.33510000E+04,
                3.12340000E+00,  1.20590000E-03, -4.20430000E-07,  6.54410000E-11, -3.76140000E-15, -1.33510000E+04,
                4.28580000E+00,  1.69940000E-04, -4.69800000E-08,  4.88050000E-12, -1.41070000E-16, -1.33510000E+04,
                2.82420000E+00, -6.30390000E-05,  4.14460000E-08, -2.12770000E-12,  3.10490000E-17, -1.33510000E+04,
               -7.28620000E+00,  1.80800000E-03, -8.78970000E-08,  1.83160000E-12, -1.42390000E-17, -1.33510000E+04,

                3.13320000E+00,  2.73250000E-03,  1.03530000E-05, -1.74770000E-08,  7.80100000E-12, -4.73900000E+04,
                4.89430000E+00,  2.34790000E-03, -8.34670000E-07,  1.31910000E-10, -7.68540000E-15, -4.73900000E+04,
                7.36680000E+00,  5.38290000E-05, -1.36000000E-08,  1.49710000E-12, -4.39190000E-17, -4.73900000E+04,
                7.16150000E+00, -5.82590000E-05,  1.58080000E-08, -7.04680000E-13,  9.43500000E-18, -4.73900000E+04,
                4.06620000E+00,  4.90030000E-04, -2.06870000E-08,  3.76430000E-13, -2.58570000E-18, -4.73900000E+04

            };

            int ns = 16, nIndex = 0, nStart = 0;
            for (int s = 0; s < ns; ++ s)
            {
                nStart = s * 30;
                for (int i = 0; i < 5; ++ i)
                {
                    nIndex = nStart + i * 6;
                    for (int j = 0; j < 6; ++ j)
                    {
                        fittingData[s][i][j] = curveFitData[nIndex + j];
                    }
                }
            }
        }

        //============================================================
        SpeciesData::SpeciesData()
        {
            nSpecies = 0;
            speciesName = NULL;
            nSpeciesType = NULL;
            nSpeciesCharge = NULL;
            speciesMass = NULL;
            nSpeciesG0 = NULL;
            nSpeciesG1 = NULL;
            speciesTtr = NULL;
            speciesTvs = NULL;
            speciesTes = NULL;
            speciesSigma = NULL;
            speciesEnthalpy = NULL;
        }

        SpeciesData::~SpeciesData()
        {
            DeAllocateMemories();
        }

        void SpeciesData::AllocateMemories()
        {
            if (this->nSpecies <= 0)
                return;
            speciesName = new string[nSpecies];
            nSpeciesType = new int[nSpecies];
            nSpeciesCharge = new int[nSpecies];
            nSpeciesG0 = new int[nSpecies];
            nSpeciesG1 = new int[nSpecies];
            speciesMass = new RDouble[nSpecies];
            speciesTtr = new RDouble[nSpecies];
            speciesTvs = new RDouble[nSpecies];
            speciesTes = new RDouble[nSpecies];
            speciesSigma = new RDouble[nSpecies];
            speciesEnthalpy = new RDouble[nSpecies];
        }

        void SpeciesData::DeAllocateMemories()
        {
            if (this->nSpecies <= 0)
                return;
            delete []speciesName;
            delete []nSpeciesType;
            delete []nSpeciesCharge;
            delete []nSpeciesG0;
            delete []nSpeciesG1;
            delete []speciesMass;
            delete []speciesTtr;
            delete []speciesTvs;
            delete []speciesTes;
            delete []speciesEnthalpy;
            delete []speciesSigma;

            nSpecies = 0;
            speciesName = NULL;
            nSpeciesType = NULL;
            nSpeciesCharge = NULL;
            speciesMass = NULL;
            nSpeciesG0 = NULL;
            nSpeciesG1 = NULL;
            speciesTtr = NULL;
            speciesTvs = NULL;
            speciesTes = NULL;
            speciesEnthalpy = NULL;
            speciesSigma = NULL;
        }

        void SpeciesData::SetSpeciesNumber(int ns)
        {
            if (ns <= 0 || ns == this->nSpecies)
                return;
            DeAllocateMemories();

            this->nSpecies = ns;
            AllocateMemories();
        }

        void SpeciesData::CreatePhysicalChemicalData(int ns/* = 11*/)
        {
            this->SetSpeciesNumber(ns);

            string *nameArray = new string[11];
            int **dataArray1 = NewPointer2 <int> (4, 11);
            RDouble **dataArray2 = NewPointer2 <RDouble> (6, 11);
            // To obtain the original data.
            this->GenerateFullPhysicalChemicalData(nameArray, dataArray1, dataArray2);

            //allocate the data according to the specific number of species.
            if (ns == 5 || ns == 7)
            {
                for (int s = 0; s < 4; ++ s)
                {
                    this->speciesName[s] = nameArray[s];

                    this->nSpeciesType[s] = dataArray1[0][s];
                    this->nSpeciesCharge[s] = dataArray1[1][s];
                    this->nSpeciesG0[s] = dataArray1[2][s];
                    this->nSpeciesG1[s] = dataArray1[3][s];

                    this->speciesMass[s] = dataArray2[0][s];
                    this->speciesTtr[s] = dataArray2[1][s];
                    this->speciesTvs[s] = dataArray2[2][s];
                    this->speciesTes[s] = dataArray2[3][s];
                    this->speciesSigma[s] = dataArray2[4][s];
                    this->speciesEnthalpy[s] = dataArray2[5][s];
                }

                if (ns == 5)
                {
                    this->speciesName[4] = nameArray[9];

                    this->nSpeciesType[4] = dataArray1[0][9];
                    this->nSpeciesCharge[4] = dataArray1[1][9];
                    this->nSpeciesG0[4] = dataArray1[2][9];
                    this->nSpeciesG1[4] = dataArray1[3][9];

                    this->speciesMass[4] = dataArray2[0][9];
                    this->speciesTtr[4] = dataArray2[1][9];
                    this->speciesTvs[4] = dataArray2[2][9];
                    this->speciesTes[4] = dataArray2[3][9];
                    this->speciesSigma[4] = dataArray2[4][9];
                    this->speciesEnthalpy[4] = dataArray2[5][9];
                }
                else
                {
                    this->speciesName[4] = nameArray[6];

                    this->nSpeciesType[4] = dataArray1[0][6];
                    this->nSpeciesCharge[4] = dataArray1[1][6];
                    this->nSpeciesG0[4] = dataArray1[2][6];
                    this->nSpeciesG1[4] = dataArray1[3][6];

                    this->speciesMass[4] = dataArray2[0][6];
                    this->speciesTtr[4] = dataArray2[1][6];
                    this->speciesTvs[4] = dataArray2[2][6];
                    this->speciesTes[4] = dataArray2[3][6];
                    this->speciesSigma[4] = dataArray2[4][6];
                    this->speciesEnthalpy[4] = dataArray2[5][6];
                    int n1 = 0, n2 = 0;
                    for (int s = 0; s < 2; ++ s)
                    {
                        n1 = s + 5;
                        n2 = s + 9;
                        this->speciesName[n1] = nameArray[n2];

                        this->nSpeciesType[n1] = dataArray1[0][n2];
                        this->nSpeciesCharge[n1] = dataArray1[1][n2];
                        this->nSpeciesG0[n1] = dataArray1[2][n2];
                        this->nSpeciesG1[n1] = dataArray1[3][n2];

                        this->speciesMass[n1] = dataArray2[0][n2];
                        this->speciesTtr[n1] = dataArray2[1][n2];
                        this->speciesTvs[n1] = dataArray2[2][n2];
                        this->speciesTes[n1] = dataArray2[3][n2];
                        this->speciesSigma[n1] = dataArray2[4][n2];
                        this->speciesEnthalpy[n1] = dataArray2[5][n2];
                    }
                }
            }
            else if (ns == 11)
            {
                for (int s = 0; s < 11; ++ s)
                {
                    this->speciesName[s] = nameArray[s];

                    this->nSpeciesType[s] = dataArray1[0][s];
                    this->nSpeciesCharge[s] = dataArray1[1][s];
                    this->nSpeciesG0[s] = dataArray1[2][s];
                    this->nSpeciesG1[s] = dataArray1[3][s];

                    this->speciesMass[s] = dataArray2[0][s];
                    this->speciesTtr[s] = dataArray2[1][s];
                    this->speciesTvs[s] = dataArray2[2][s];
                    this->speciesTes[s] = dataArray2[3][s];
                    this->speciesSigma[s] = dataArray2[4][s];
                    this->speciesEnthalpy[s] = dataArray2[5][s];
                }
            }

            delete []nameArray;
            DelPointer2(dataArray1);
            DelPointer2(dataArray2);
        }

        void SpeciesData::GenerateFullPhysicalChemicalData(string *nameArray, int **dataArray1, RDouble **dataArray2)
        {
            // The name of species.
            nameArray[0] = "O";
            nameArray[1] = "O2";
            nameArray[2] = "NO";
            nameArray[3] = "N";
            nameArray[4] = "O+";
            nameArray[5] = "O2+";
            nameArray[6] = "NO+";
            nameArray[7] = "N+";
            nameArray[8] = "N2+";
            nameArray[9] = "N2";
            nameArray[10] = "e-";
            // The type and charge of species.
            int data1[44] = {
                1,  2,  2,  1,  1,  2,  2,  1,  2,  2,  0,
                0,  0,  0,  0,  1,  1,  1,  1,  1,  0, -1,
                9,  3,  4,  9,  4,  4,  1,  9,  2,  1,  0,
                5,  2,  8, 10, 10, 10,  3,  5,  4,  3,  0
            };
            // The other parameter.
            RDouble data2[66] = {
                15.9994, 31.9988, 30.00614, 14.00674, 15.99885, 31.99825, 30.00564, 14.00619, 28.01293, 28.01348, 5.486e-4,
                216.0, 37.0, 106.7, 260.0, 216.0, 37.0, 106.7, 260.0, 106.7, 106.7, 106.7,
                0.0, 2239.0, 2699.0, 0.0, 0.0, 2652.0, 3373.0, 0.0, 3129.0, 3353.0, 0.0,
                22890.0, 11390.0, 55874.0, 27670.0, 38610.0, 47460.0, 75140.0, 22052.0, 13200.0, 72225.0, 0.0,
                3.1177e-19, 7.0149e-20, 1.6178e-19, 3.2471e-19, 1.5151e-15, 5.1636e-14, 6.2090e-15, 1.0457e-15, 6.2090e-15, 1.8273e-19, 6.2090e-15,
                1.556308E+07, 0.0, 3.012717E+06, 3.376946E+07, 9.813204E+07, 3.656450E+07, 3.309383E+07, 1.342264E+08, 5.461763E+07, 0.0, 0.0
            };

            int nIndex = 0;
            for (int i = 0; i < 4; ++ i)
            {
                nIndex = i * 11;
                for (int s = 0; s < 11; ++ s)
                {
                    dataArray1[i][s] = data1[nIndex + s];
                }
            }
            for (int i = 0; i < 6; ++ i)
            {
                nIndex = i * 11;
                for (int s = 0; s < 11; ++ s)
                {
                    dataArray2[i][s] = data2[nIndex + s];
                }
            }
        }

        int SpeciesData::GetSpeciesIndex(string speciesName)
        {
            int nIndex = 0;
            if (speciesName == "O")
            {
                nIndex = 0;
            }
            else if (speciesName == "O2")
            {
                nIndex = 1;
            }
            else if (speciesName == "NO")
            {
                nIndex = 2;
            }
            else if (speciesName == "N")
            {
                nIndex = 3;
            }
            else if (speciesName == "O+")
            {
                nIndex = 4;
            }
            else if (speciesName == "O2+")
            {
                nIndex = 5;
            }
            else if (speciesName == "NO+")
            {
                nIndex = 6;
            }
            else if (speciesName == "N+")
            {
                nIndex = 7;
            }
            else if (speciesName == "N2+")
            {
                nIndex = 8;
            }
            else if (speciesName == "N2")
            {
                nIndex = 9;
            }
            else if (speciesName == "e-")
            {
                nIndex = 10;
            }
            return nIndex;
        }

        string SpeciesData::GetSpeciesName(int nIndex)
        {
            string speciesName = "";
            if (nIndex == 0)
            {
                speciesName = "O";
            }
            else if (nIndex == 1)
            {
                speciesName = "O2";
            }
            else if (nIndex == 2)
            {
                speciesName = "NO";
            }
            else if (nIndex == 3)
            {
                speciesName = "N";
            }
            else if (nIndex == 4)
            {
                speciesName = "O+";
            }
            else if (nIndex == 5)
            {
                speciesName = "O2+";
            }
            else if (nIndex == 6)
            {
                speciesName = "NO+";
            }
            else if (nIndex == 7)
            {
                speciesName = "N+";
            }
            else if (nIndex == 8)
            {
                speciesName = "N2+";
            }
            else if (nIndex == 9)
            {
                speciesName = "N2";
            }
            else if (nIndex == 10)
            {
                speciesName = "e-";
            }
            return speciesName;
        }

        //=======================================
        TransportData::TransportData()
        {
            nSpecies = 0; //number of species;
            viscosityFitData = NULL; //the curve fitting data of viscosity.
            collisionAreaFitData = NULL; // the curve fitting data of the average collision area.
        }

        TransportData::~TransportData()
        {
            DeAllocateMemories();
        }

        void TransportData::AllocateMemories()
        {
            if (this->nSpecies <= 0)
                return;
            viscosityFitData = NewPointer2 <RDouble> (this->nSpecies, 5);
            collisionAreaFitData = NewPointer2 <RDouble> (this->nSpecies, 4);
        }

        void TransportData::DeAllocateMemories()
        {
            if (this->nSpecies <= 0)
                return;
            DelPointer2(viscosityFitData);
            DelPointer2(collisionAreaFitData);
            nSpecies = 0;
            viscosityFitData = NULL;
            collisionAreaFitData = NULL;
        }

        void TransportData::SetSpeciesNumber(int ns)
        {
            if (ns <= 0 || ns == this->nSpecies)
                return;
            DeAllocateMemories();

            this->nSpecies = ns;
            AllocateMemories();
        }

        void TransportData::CreateCurveFittingData(int nGasModel/* = 0*/, int ns/* = 11*/)
        {
            this->SetSpeciesNumber(ns);

            RDouble **viscosityData, **collisionData;
            if (nGasModel == 0)
            {
                viscosityData = NewPointer2 <RDouble> (11, 5);
                collisionData = NewPointer2 <RDouble> (11, 4);
            }
            else
            {
                viscosityData = NewPointer2 <RDouble> (12, 5);
                collisionData = NewPointer2 <RDouble> (12, 4);
            }

            // To generate the initial data of curve fitting method.
            this->GenerateFullCurveFittingData(viscosityData, collisionData, nGasModel);

            // Allocate the curve fitting data according to the number of species.
            if (nGasModel == 0) //The Earth Atmosphere.
            {
                if (ns == 5 || ns == 7)
                {
                    for (int s = 0; s < 4; ++ s)
                    {
                        for (int i = 0; i < 5; ++ i)
                        {
                            this->viscosityFitData[s][i] = viscosityData[s][i];
                        }
                        for (int i = 0; i < 4; ++ i)
                        {
                            this->collisionAreaFitData[s][i] = collisionData[s][i];
                        }
                    }
                    if (ns == 5)
                    {
                        for (int i = 0; i < 5; ++ i)
                        {
                            this->viscosityFitData[4][i] = viscosityData[9][i];
                        }
                        for (int i = 0; i < 4; ++ i)
                        {
                            this->collisionAreaFitData[4][i] = collisionData[9][i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < 5; ++ i)
                        {
                            this->viscosityFitData[4][i] = viscosityData[6][i];
                            this->viscosityFitData[5][i] = viscosityData[9][i];
                            this->viscosityFitData[6][i] = viscosityData[10][i];
                        }
                        for (int i = 0; i < 4; ++ i)
                        {
                            this->collisionAreaFitData[4][i] = collisionData[6][i];
                            this->collisionAreaFitData[5][i] = collisionData[9][i];
                            this->collisionAreaFitData[6][i] = collisionData[10][i];
                        }
                    }
                }
                else if (ns == 11)
                {
                    for (int s = 0; s < 11; ++ s)
                    {
                        for (int i = 0; i < 5; ++ i)
                        {
                            this->viscosityFitData[s][i] = viscosityData[s][i];
                        }
                        for (int i = 0; i < 4; ++ i)
                        {
                            this->collisionAreaFitData[s][i] = collisionData[s][i];
                        }
                    }
                }
            }
            else //The Mars Atmosphere.
            {
                if (ns == 8 || ns == 12)
                {
                    for (int s = 0; s < ns; ++ s)
                    {
                        for (int i = 0; i < 5; ++ i)
                        {
                            this->viscosityFitData[s][i] = viscosityData[s][i];
                        }
                        for (int i = 0; i < 4; ++ i)
                        {
                            this->collisionAreaFitData[s][i] = collisionData[s][i];
                        }
                    }
                }
            }

            DelPointer2(viscosityData);
            DelPointer2(collisionData);
        }

        void TransportData::GenerateFullCurveFittingData(RDouble **viscosityData, RDouble **collisionData, int nGasModel/* = 0*/)
        {
            int nIndex = 0;
            int speciesNumber = 11;
            if (nGasModel == 0) // The Earth atmosphere.
            {
                speciesNumber = 11;
                // the order of species is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
                RDouble curveFitData1[55] = {
                    0.0,      0.0,    0.0203144,  0.4294404, -11.6031403,
                    0.0,      0.0,    0.044929,  -0.0826158, -9.2019475,
                    0.0,      0.0,    0.0436378, -0.0335511, -9.5767430,
                    0.0,      0.0,    0.0115572,  0.6031679, -12.432750,
                    0.0912,  -3.3154,  45.1290,  -270.4211,  586.2903,
                    0.0908,  -3.3020,  44.9511,  -269.3877,  584.4130,
                    0.0913,  -3.3178,  45.1426,  -270.3948,  586.3300,
                    0.0895,  -3.2573,  44.3511,  -265.8276,  576.1313,
                    0.0897,  -3.2618,  44.4079,  -266.1462,  577.1449,
                    0.0,      0.0,     0.0268142,  0.3177838,  -11.3155513,
                    0.0899,  -3.2731,  44.5782,  -267.2522,  574.4149,
                };
                RDouble curveFitData2[44] = {
                    0.0000, -0.0207,  0.0780,  3.5658,
                    0.0000, -0.0485,  0.6475,  1.2607,
                    0.0000, -0.0453,  0.5624,  1.7669,
                    0.0000, -0.0118, -0.0960,  4.3252,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000, -0.0203,  0.0683,  4.0900,
                    0.0000,  0.0000, -2.0000,  24.3061
                };

                for (int s = 0; s < speciesNumber; ++ s)
                {
                    nIndex = s * 5;
                    for (int i = 0; i < 5; ++ i)
                    {
                        viscosityData[s][i] = curveFitData1[nIndex + i];
                    }

                    nIndex = s * 4;
                    for (int i = 0; i < 4; ++ i)
                    {
                        collisionData[s][i] = curveFitData2[nIndex + i];
                    }
                }

            }
            else // The Mars atmosphere.
            {
                speciesNumber = 12;
                // the order of species is O, O2, NO, N, N2, C, CO, CO2, C2, CN, NO+, e-.
                RDouble curveFitData1[60] = {
                    0.0000000,   0.0000000,   0.0203144,   0.4294404,   -11.6031400,
                    0.0000000,   0.0000000,   0.0449290,  -0.0826158,   -9.2019475,
                    0.0000000,   0.0000000,   0.0436378,  -0.0335511,   -9.5767430,
                    0.0000000,   0.0000000,   0.0115572,   0.6031679,   -12.4327500,
                    0.0000000,   0.0000000,   0.0268142,   0.3177838,   -11.3155510,
                    0.0000000,   0.0000000,  -0.0061100,   0.7385000,   -12.7100000,
                    0.0000000,   0.0000000,  -0.0324200,   1.1430000,   -14.1700000,
                    0.0000000,   0.0000000,  -0.0577200,   1.5430000,   -15.7600000,
                    0.0000000,   0.0000000,  -0.0272900,   1.0640000,   -14.0300000,
                    0.0000000,   0.0000000,  -0.0257200,   1.0400000,   -13.8600000,
                    0.0000000,   0.0000000,   0.3020141,  -3.5039791,   -3.7355157,
                    0.0000000,   0.0000000,   0.0000000,   0.0000000,   -20.0000000
                };
                RDouble curveFitData2[48] = {
                    0.0000, -0.0207,  0.0780,  3.5658,
                    0.0000, -0.0485,  0.6475,  1.2607,
                    0.0000, -0.0453,  0.5624,  1.7669,
                    0.0000, -0.0118, -0.0960,  4.3252,
                    0.0000, -0.0203,  0.0683,  4.0900,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3602,
                    0.0000,  0.0000, -2.0000,  24.3061
                };

                for (int s = 0; s < speciesNumber; ++ s)
                {
                    nIndex = s * 5;
                    for (int i = 0; i < 5; ++ i)
                    {
                        viscosityData[s][i] = curveFitData1[nIndex + i];
                    }

                    nIndex = s * 4;
                    for (int i = 0; i < 4; ++ i)
                    {
                        collisionData[s][i] = curveFitData2[nIndex + i];
                    }
                }

            }//End of else expression.

        }

        //============================================================
        MarsSpeciesData::MarsSpeciesData()
        {
            nSpecies = 0;
            speciesName = NULL;
            nSpeciesType = NULL;
            nSpeciesCharge = NULL;
            speciesMass = NULL;
            nSpeciesG0 = NULL;
            nSpeciesG1 = NULL;
            speciesTtr = NULL;
            speciesTvs = NULL;
            speciesTes = NULL;
            speciesSigma = NULL;
            speciesEnthalpy = NULL;
        }

        MarsSpeciesData::~MarsSpeciesData()
        {
            DeAllocateMemories();
        }

        void MarsSpeciesData::AllocateMemories()
        {
            if (this->nSpecies <= 0)
                return;
            speciesName = new string[nSpecies];
            nSpeciesType = new int[nSpecies];
            nSpeciesCharge = new int[nSpecies];
            nSpeciesG0 = new int[nSpecies];
            nSpeciesG1 = new int[nSpecies];
            speciesMass = new RDouble[nSpecies];
            speciesTtr = new RDouble[nSpecies];
            speciesTvs = new RDouble[nSpecies];
            speciesTes = new RDouble[nSpecies];
            speciesSigma = new RDouble[nSpecies];
            speciesEnthalpy = new RDouble[nSpecies];
        }

        void MarsSpeciesData::DeAllocateMemories()
        {
            if (this->nSpecies <= 0)
                return;
            delete []speciesName;
            delete []nSpeciesType;
            delete []nSpeciesCharge;
            delete []nSpeciesG0;
            delete []nSpeciesG1;
            delete []speciesMass;
            delete []speciesTtr;
            delete []speciesTvs;
            delete []speciesTes;
            delete []speciesEnthalpy;
            delete []speciesSigma;

            nSpecies = 0;
            speciesName = NULL;
            nSpeciesType = NULL;
            nSpeciesCharge = NULL;
            speciesMass = NULL;
            nSpeciesG0 = NULL;
            nSpeciesG1 = NULL;
            speciesTtr = NULL;
            speciesTvs = NULL;
            speciesTes = NULL;
            speciesEnthalpy = NULL;
            speciesSigma = NULL;
        }

        void MarsSpeciesData::SetSpeciesNumber(int ns)
        {
            if (ns <= 0 || ns == this->nSpecies)
                return;
            DeAllocateMemories();

            this->nSpecies = ns;
            AllocateMemories();
        }

        void MarsSpeciesData::CreatePhysicalChemicalData(int ns/* = 12*/)
        {
            if (ns != 8 && ns != 12)
            {
                return;
            }

            this->SetSpeciesNumber(ns);

            string *nameArray = new string[12];
            int **dataArray1 = NewPointer2 <int> (4, 12);
            RDouble **dataArray2 = NewPointer2 <RDouble> (6, 12);
            // To obtain the original data.
            this->GenerateFullPhysicalChemicalData(nameArray, dataArray1, dataArray2);

            //allocate the data according to the specific number of species.
            for (int s = 0; s < ns; ++ s)
            {
                this->speciesName[s] = nameArray[s];

                this->nSpeciesType[s] = dataArray1[0][s];
                this->nSpeciesCharge[s] = dataArray1[1][s];
                this->nSpeciesG0[s] = dataArray1[2][s];
                this->nSpeciesG1[s] = dataArray1[3][s];

                this->speciesMass[s] = dataArray2[0][s];
                this->speciesTtr[s] = dataArray2[1][s];
                this->speciesTvs[s] = dataArray2[2][s];
                this->speciesTes[s] = dataArray2[3][s];
                this->speciesSigma[s] = dataArray2[4][s];
                this->speciesEnthalpy[s] = dataArray2[5][s];
            }

            delete []nameArray;
            DelPointer2(dataArray1);
            DelPointer2(dataArray2);
        }

        void MarsSpeciesData::GenerateFullPhysicalChemicalData(string *nameArray, int **dataArray1, RDouble **dataArray2)
        {
            // The name of species.
            nameArray[0] = "O";
            nameArray[1] = "O2";
            nameArray[2] = "NO";
            nameArray[3] = "N";
            nameArray[4] = "N2";
            nameArray[5] = "C";
            nameArray[6] = "CO";
            nameArray[7] = "CO2";
            nameArray[8] = "C2";
            nameArray[9] = "CN";
            nameArray[10] = "NO+";
            nameArray[11] = "e-";
            // The type and charge of species.
            int data1[48] = {
                1,  2,  2,  1,  2,  1,  2,  2,  2,  2,  2,  0,
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,
                9,  3,  4,  9,  1,  9,  1,  1,  6,  2,  1,  0,
                5,  2,  8, 10,  3,  5,  6,  2,  6,  4,  3,  0
            };

            // The other parameter.
            RDouble data2[72] = {
                15.9994, 31.9988, 30.00614, 14.00674, 28.01293, 12.011, 28.01, 44.01, 24.022, 26.01774, 30.00564, 5.486e-4,
                216.0, 37.0, 106.7, 260.0, 106.7, 216.0, 216.0, 216.0, 216.0, 216.0, 106.7, 106.7,
                0.0, 2239.0, 2817.0, 0.0, 3395, 0.0, 3074, 945.0, 2362.0, 2976, 3373.0, 0.0,
                22890.0, 11390.0, 55874.0, 27670.0, 72225.0, 14670.0, 70060.0, 66200.0, 27800.0, 13300.0, 75140.0, 0.0,
                3.1177e-19, 7.0149e-19, 1.6178e-19, 3.2471e-19, 1.8273e-19, 3.1177e-19, 3.1177e-19, 3.1177e-19, 3.1177e-19, 3.1177e-19, 6.2090e-19, 6.2090e-24,
                1.556308E+07, 0.0, 3.012717E+06, 3.376946E+07, 0.0, 5.919574E+07, -4.069975E+06, -8.929789E+06, 3.421863E+07, 1.679623E+07, 3.309383E+07, 0.0
            };

            int nIndex = 0;
            int speciesNumber = 12;
            for (int i = 0; i < 4; ++ i)
            {
                nIndex = i * speciesNumber;
                for (int s = 0; s < speciesNumber; ++ s)
                {
                    dataArray1[i][s] = data1[nIndex + s];
                }
            }
            for (int i = 0; i < 6; ++ i)
            {
                nIndex = i * speciesNumber;
                for (int s = 0; s < speciesNumber; ++ s)
                {
                    dataArray2[i][s] = data2[nIndex + s];
                }
            }
        }

        int MarsSpeciesData::GetSpeciesIndex(string speciesName)
        {
            int nIndex = 0;
            if (speciesName == "O")
            {
                nIndex = 0;
            }
            else if (speciesName == "O2")
            {
                nIndex = 1;
            }
            else if (speciesName == "NO")
            {
                nIndex = 2;
            }
            else if (speciesName == "N")
            {
                nIndex = 3;
            }
            else if (speciesName == "N2")
            {
                nIndex = 4;
            }
            else if (speciesName == "C")
            {
                nIndex = 5;
            }
            else if (speciesName == "CO")
            {
                nIndex = 6;
            }
            else if (speciesName == "CO2")
            {
                nIndex = 7;
            }
            else if (speciesName == "C2")
            {
                nIndex = 8;
            }
            else if (speciesName == "CN")
            {
                nIndex = 9;
            }
            else if (speciesName == "NO+")
            {
                nIndex = 10;
            }
            else if (speciesName == "e-")
            {
                nIndex = 11;
            }
            return nIndex;
        }

        string MarsSpeciesData::GetSpeciesName(int nIndex)
        {
            string speciesName = "";
            if (nIndex == 0)
            {
                speciesName = "O";
            }
            else if (nIndex == 1)
            {
                speciesName = "O2";
            }
            else if (nIndex == 2)
            {
                speciesName = "NO";
            }
            else if (nIndex == 3)
            {
                speciesName = "N";
            }
            else if (nIndex == 4)
            {
                speciesName = "N2";
            }
            else if (nIndex == 5)
            {
                speciesName = "C";
            }
            else if (nIndex == 6)
            {
                speciesName = "CO";
            }
            else if (nIndex == 7)
            {
                speciesName = "CO2";
            }
            else if (nIndex == 8)
            {
                speciesName = "C2";
            }
            else if (nIndex == 9)
            {
                speciesName = "CN";
            }
            else if (nIndex == 10)
            {
                speciesName = "NO+";
            }
            else if (nIndex == 11)
            {
                speciesName = "e-";
            }
            return speciesName;
        }

    }// end of GAS_SPACE

}// end of PHSPACE
