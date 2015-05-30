#ifndef __LASM_TerminatorChemistryTest__
#define __LASM_TerminatorChemistryTest__

#ifdef LASM_IN_SPHERE

#include "DeformationTest.h"

class TerminatorChemistryTest : public DeformationTest {
protected:
    SpaceCoord k1Center;
    const double Xt0 = 4e-6;
public:
    TerminatorChemistryTest();
    virtual ~TerminatorChemistryTest();

    virtual void
    init(const ConfigManager &configManager, AdvectionManager &advactionManager);

    virtual void
    setInitialCondition(AdvectionManager &advactionManager);

    virtual void
    advancePhysics(const TimeLevelIndex<2> &timeIdx,
                   AdvectionManager &advectionManager);

    virtual void
    output(const TimeLevelIndex<2> &timeIdx,
           const AdvectionManager &advactionManager);
protected:
    void
    calcChemicalReactionRate(const SpaceCoord &x, double &k1, double &k2);

    void
    calcTendency(const SpaceCoord &x, double X, double X2,
                 double &dX, double &dX2);
}; // TerminatorChemistryTest

#endif // LASM_IN_SPHERE

#endif // __LASM_TerminatorChemistryTest__
