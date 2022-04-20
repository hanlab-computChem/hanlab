//
//  main.cpp
//  Adodemo
//
//  Created by Xuan Tang on 2020/2/24.
//  Copyright © 2020 tangx. All rights reserved.
//

#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Angle.h"
#include "tools/Pbc.h"
#include "core/PlumedMain.h"

#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;

#define PI acos(-1)

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR NAME
/*
*/
//+ENDPLUMEDOC

/**We begin by declaring a class for your colvar.  This class inherits everything from the Colvar class.
       This ensures it has a label, a place to store its value, places to the store the values of the derivatives
       and that it can access the various atoms it will employ.
*/

class ADO : public Colvar {
public:
    static void registerKeywords( Keywords& keys );
    explicit ADO(const ActionOptions&);
    virtual void calculate();
protected:
    //bool pbc;
    double R_0, D_0, SigOrien, SigOneD1, reforient, SigOrien2, reforient2, refangle1, SigOneD2, refangle2, dmax;
private:
    unsigned int AtomNum;
    int size1, size2;
};

// The following command inserts your new colvar into plumed by inserting calls to your new routines into the parts of plumed where they are required.  This macro takes two arguments:
// The first is the name of your ColvarClass and the second is the keyword for your CV (the first word in the input line for your CV).

PLUMED_REGISTER_ACTION(ADO,"ADO")

// The following routine creates the documentation for the keyowrds used by your CV
void ADO::registerKeywords( Keywords& keys )
{
    Colvar::registerKeywords(keys);
    keys.add("compulsory","NN","6","The n parameter of the switching function");
    keys.add("compulsory","MM","12","The m parameter of the switching function; 0 implies 2*NN");
    keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
    keys.add("compulsory","R_0","The r_0 parameter of the switching function");
    keys.add("compulsory","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. ");
    
    keys.add("compulsory", "REFO", "The reference angle of orientaion");
    keys.add("compulsory", "REFA1", "The reference angle of one-dimensional arrangement");
    
    keys.add("compulsory", "SIGORIEN", "0.4", "The sigma value of orientaion");
    keys.add("compulsory", "SIGONED1", "0.4", "The sigma value of one-dimensional arrangement");
    keys.add("atoms","GROUPA","First list of atoms");
    keys.add("atoms","GROUPB","Second list of atoms");
    keys.add("optional", "REFO2", "The second angle of orientaion");
    keys.add("optional", "SIGORIEN2", "The second sigma value of orientaion");
    keys.add("optional", "REFA2", "The second reference angle of one-dimensional arrangement");
    keys.add("optional", "SIGONED2", "The second sigma value of one-dimensional arrangement");
    keys.add("optional","D_MAX","The maximum of having contact");
}

// We now write the actual readin (constructor) and calculations routines for the colvar

ADO::ADO(const ActionOptions&ao):
// This line sets up various things in the plumed core which colvars rely on.
    PLUMED_COLVAR_INIT(ao)

    //pbc(true)
{
    parse("R_0",R_0);
    parse("D_0",D_0);
    parse("SIGORIEN", SigOrien);
    parse("SIGONED1", SigOneD1);
    parse("SIGORIEN2", SigOrien2);
    parse("SIGONED2", SigOneD2);
    parse("D_MAX", dmax);
    parse("REFO", reforient);
    parse("REFA1", refangle1);
    parse("REFA2", refangle2);
    parse("REFO2", reforient2);
    vector<AtomNumber> atoms;
    vector<AtomNumber> atoms1;
    vector<AtomNumber> atoms2; // You almost always have atoms
    parseAtomList("GROUPA", atoms1);
    parseAtomList("GROUPB", atoms2);
    AtomNum = atoms1.size();
    
    if(AtomNum < 1)
    {
        error("Number of specified atoms should be larger than 0");
    }
    /*
    bool nopbc=!pbc;
    parseFlag("NOPBC",nopbc);
    pbc=!nopbc;
    */
    //calculate atoms
    size1 = atoms1.size();
    size2 = atoms2.size();
    if(atoms2.size()>0)
    {
        for(int i = 0; i < size1; i++)
        {
            AtomNumber a1;
            a1.setIndex(atoms1[i].index());
            atoms.push_back(a1);
        }
        
        for(int i = 0; i < size2; i++)
        {
            AtomNumber a2;
            a2.setIndex(atoms2[i].index());
            atoms.push_back(a2);
        }
    }
    else
    {
        for(int i = 0; i < size1; i++)
        {
            AtomNumber a1;
            a1.setIndex(atoms1[i].index());
            atoms.push_back(a1);
        }
    }
    checkRead(); //This command checks that everything on the input line has been read properly
    
    //The following two lines inform the plumed core that we require space to store the value of the CV and that the CV will act on a particular list of atoms.
    addValueWithDerivatives(); setNotPeriodic();
    
    requestAtoms(atoms);
}

void ADO::calculate()
{
    //if(pbc) makeWhole();
    vector<Vector> pos = getPositions();
    int numTotal = pos.size();
    vector<Vector> derivatives;
    double cv_val = 0; 
    Tensor virial;
    //vector<Vector> der_val(numTotal);
    if(size2 > 0)
    {
		vector<Vector> der1(size1/2); //c1
        vector<Vector> der2(size1/2); //o1
        vector<Vector> der3(size2/2); //c2
        vector<Vector> der4(size2/2); //o2
        //These are the things you must calculate for any cv
        for(int i = 0; i < size1/2; i++)
        {
        
            for(int j = 0; j < size2/2; j++)
            {
                if(i != j)
                {
                    Vector dist = pbcDistance(pos[i],pos[j+size1]); //C1-C2
                    double distlen = dist.modulo();

					Vector dis1 = pbcDistance(pos[i], pos[i+size1/2]); // C1-O1
                    Vector dis2 = pbcDistance(pos[j+size1], pos[j+size1+size2/2]); // C2-O2
                    
                    Vector dOrien1, dOrien2;
                    Vector dAngle1, dAngle2;
                    PLMD::Angle a;
                    
                    double Orien, Angle;
                    Orien = a.compute(dis1, dis2, dOrien1, dOrien2);
                    Angle = a.compute(dis1, dist, dAngle1, dAngle2);

                    if(distlen > dmax)
                    {
                    	cv_val += 0.;
                    	der1[i] += 0*(-dOrien1);
	                    der2[i] += 0*dOrien1;
	                    der3[j] += 0*(-dOrien2);
	                    der4[j] += 0*dOrien2;
	                    
	                    Vector ddfunc1, ddfunc2, ddfunc3, ddfunc4, ddfunc5;
	                    ddfunc1 = 0*dOrien1;
	                    Tensor vv1(ddfunc1,dis1);
	                    virial -= vv1;
	                    
	                    ddfunc2 = 0*dOrien2;
	                    Tensor vv2(ddfunc2,dis2);
	                    virial -= vv2;
	                    
	                    ddfunc3 = 0*dAngle1;
	                    Tensor vv3(ddfunc3,dis1);
	                    virial -= vv3;
	                    
	                    ddfunc4 = 0*dAngle2;
	                    Tensor vv4(ddfunc4,dist);
	                    virial -= vv4;
	                    
	                    ddfunc5 = 0*dist;
	                    Tensor vv5(ddfunc5,dist);
	                    virial -= vv5;
                    }
                    else
                    {
	                    double scoreOrien;
	                    double scoreAngle;
	                    double scoreDis;
	                    
	                    double pre_Orien;

			            if(SigOrien2)
			            {
			            	double rorient1 = (cos(Orien)-reforient)/SigOrien;
			            	double rorient2 = (cos(Orien)-reforient2)/SigOrien2;
			            	scoreOrien = exp(-0.5*rorient1*rorient1) + exp(-0.5*rorient2*rorient2);
			            	pre_Orien = -exp(-0.5*rorient1*rorient1)/SigOrien*rorient1*(-sin(Orien)) - exp(-0.5*rorient2*rorient2)/SigOrien2*rorient2*(-sin(Orien));
			            }
			            else
			            {
			            	double rorient = (cos(Orien)-reforient)/SigOrien;
			            	scoreOrien = exp(-0.5*rorient*rorient);
			            	pre_Orien = -scoreOrien/SigOrien*rorient*(-sin(Orien));
			            }

	                    double pre_Angle;

	                    if(SigOneD2)
	                    {
	                        double rangle1 = (cos(Angle)-refangle1)/SigOneD1;
	                        double rangle2 = (cos(Angle)-refangle2)/SigOneD2;
	                        scoreAngle = exp(-0.5*rangle1*rangle1) + exp(-0.5*rangle2*rangle2);
	                        pre_Angle = -exp(-0.5*rangle1*rangle1)/SigOneD1*rangle1*(-sin(Angle)) - exp(-0.5*rangle2*rangle2)/SigOneD2*rangle2*(-sin(Angle));
	                    }
	                    else
	                    {
	                        double rangle = (cos(Angle)-refangle1)/SigOneD1;                        
	                        scoreAngle = exp(-0.5*rangle*rangle);
	                        pre_Angle = -scoreAngle/SigOneD1*rangle * (-sin(Angle));
	                    }
	                    
	                    double pre_Dist;
	                    
			    double rdis = (distlen-D_0)/R_0;
	                    scoreDis = exp(-0.5*rdis*rdis);
	                    pre_Dist = -scoreDis*rdis/R_0;           

	                    //The value of the cv
	                    cv_val += scoreOrien * scoreAngle * scoreDis / (numTotal/2/2);
	                    
	                    Vector der_Dist;
	                    der_Dist = pre_Dist * dist / distlen;
	                    
	                    der1[i] += (pre_Orien*(-dOrien1)*scoreAngle*scoreDis + pre_Angle*(-dAngle1)*scoreOrien*scoreDis + pre_Angle*(-dAngle2)*scoreOrien*scoreDis + (-der_Dist*scoreOrien*scoreAngle))/(numTotal/2/2);

	                    der2[i] += (pre_Orien*dOrien1*scoreAngle*scoreDis + pre_Angle*dAngle1*scoreOrien*scoreDis)/(numTotal/2/2);
	                    
	                    der3[j] += (pre_Orien*(-dOrien2)*scoreAngle*scoreDis + pre_Angle*dAngle2*scoreOrien*scoreDis + der_Dist*scoreOrien*scoreAngle)/(numTotal/2/2);
	                    
	                    der4[j] += pre_Orien*dOrien2*scoreAngle*scoreDis/(numTotal/2/2);
	                    
	                    Vector ddfunc1, ddfunc2, ddfunc3, ddfunc4, ddfunc5;
	                    ddfunc1 = pre_Orien*dOrien1*scoreAngle*scoreDis/(numTotal/2/2);
	                    Tensor vv1(ddfunc1,dis1);
	                    virial -= vv1;
	                    
	                    ddfunc2 = pre_Orien*dOrien2*scoreAngle*scoreDis/(numTotal/2/2);
	                    Tensor vv2(ddfunc2,dis2);
	                    virial -= vv2;
	                    
	                    ddfunc3 = pre_Angle*dAngle1*scoreOrien*scoreDis/(numTotal/2/2);
	                    Tensor vv3(ddfunc3,dis1);
	                    virial -= vv3;
	                    
	                    ddfunc4 = pre_Angle*dAngle2*scoreOrien*scoreDis/(numTotal/2/2);
	                    Tensor vv4(ddfunc4,dist);
	                    virial -= vv4;
	                    
	                    ddfunc5 = der_Dist*scoreOrien*scoreAngle/(numTotal/2/2);
	                    Tensor vv5(ddfunc5,dist);
	                    virial -= vv5;
	                }
                }
            }
        }
	for(int i = 0; i < der1.size(); i++)
	{
		derivatives.push_back(der1[i]);
	}
	derivatives.insert(derivatives.end(), der2.begin(), der2.end());
	derivatives.insert(derivatives.end(), der3.begin(), der3.end());
	derivatives.insert(derivatives.end(), der4.begin(), der4.end());
    }
    else
    {
		vector<Vector> der1(size1); 
        for(int i = 0; i < numTotal/2-1; i++)
        {
            for(int j = i+1; j < numTotal/2; j++)
            {
                Vector dist = pbcDistance(pos[i], pos[j]);
                double distlen = dist.modulo();
            	Vector dis1 = pbcDistance(pos[i], pos[i+numTotal/2]);
                Vector dis2 = pbcDistance(pos[j], pos[j+numTotal/2]);
                
                Vector dOrien1, dOrien2;
                Vector dAngle1, dAngle2;
                PLMD::Angle a;
                
                double Orien, Angle;
                Orien = a.compute(dis1, dis2, dOrien1, dOrien2);
                Angle = a.compute(dis1, dist, dAngle1, dAngle2);
                if(distlen > dmax)
                {
                	cv_val += 0.;
                	der1[i] += 0*(-dOrien1);
	                der1[j] += 0*dOrien1;
	                der1[i+numTotal/2] += 0*(-dOrien2);
	                der1[j+numTotal/2] += 0*dist;
	                
	                Vector ddfunc1, ddfunc2, ddfunc3, ddfunc4, ddfunc5;
	                ddfunc1 = 0*dOrien1;
	                Tensor vv1(ddfunc1,dis1);
	                virial -= vv1;
	                
	                ddfunc2 = 0*dOrien2;
	                Tensor vv2(ddfunc2,dis2);
	                virial -= vv2;
	                
	                ddfunc3 = 0*dAngle1;
	                Tensor vv3(ddfunc3,dis1);
	                virial -= vv3;
	                
	                ddfunc4 = 0*dAngle2;
	                Tensor vv4(ddfunc4,dist);
	                virial -= vv4;
	                
	                ddfunc5 = 0*dist;
	                Tensor vv5(ddfunc5,dist);
	                virial -= vv5;
                }
                else
                {   
	                double scoreOrien;
	                double scoreAngle;
	                double scoreDis;

	                double pre_Orien;

		            if(SigOrien2)
		            {
		            	double rorient1 = (cos(Orien)-reforient)/SigOrien;
		            	double rorient2 = (cos(Orien)-reforient2)/SigOrien2;
		            	scoreOrien = exp(-0.5*rorient1*rorient1) + exp(-0.5*rorient2*rorient2);
		            	pre_Orien = -exp(-0.5*rorient1*rorient1)/SigOrien*rorient1*(-sin(Orien)) - exp(-0.5*rorient2*rorient2)/SigOrien2*rorient2*(-sin(Orien));
		            }
		            else
		            {
		            	double rorient = (cos(Orien)-reforient)/SigOrien;
		            	scoreOrien = exp(-0.5*rorient*rorient);
		            	pre_Orien = -scoreOrien/SigOrien*rorient*(-sin(Orien));
		            }

	                double pre_Angle;

	                if(SigOneD2)
	                {
	                    double rangle1 = (cos(Angle)-refangle1)/SigOneD1;
	                    double rangle2 = (cos(Angle)-refangle2)/SigOneD2;
	                    scoreAngle = exp(-0.5*rangle1*rangle1) + exp(-0.5*rangle2*rangle2);
	                    pre_Angle = -exp(-0.5*rangle1*rangle1)/SigOneD1*rangle1*(-sin(Angle)) - exp(-0.5*rangle2*rangle2)/SigOneD2*rangle2*(-sin(Angle));
	                }
	                else
	                {
	                    double rangle = (cos(Angle)-refangle1)/SigOneD1;                        
	                    scoreAngle = exp(-0.5*rangle*rangle);
	                    pre_Angle = -scoreAngle/SigOneD1*rangle * (-sin(Angle));
	                }
	                
	                double pre_Dist;
	                
	                double rdis = (distlen-D_0)/R_0;
	                scoreDis = exp(-0.5*rdis*rdis);
	                pre_Dist = -scoreDis*rdis/R_0; 

	                //The value of the cv
	                cv_val += scoreOrien * scoreAngle * scoreDis / (numTotal/2);
	                
	                Vector der_Dist;
	                der_Dist = pre_Dist * dist / distlen;
	                
	                der1[i] += (pre_Orien*(-dOrien1)*scoreAngle*scoreDis + pre_Angle*(-dAngle1)*scoreOrien*scoreDis + pre_Angle*(-dAngle2)*scoreOrien*scoreDis + (-der_Dist*scoreOrien*scoreAngle))/(numTotal/2);

	                der1[j] += (pre_Orien*dOrien1*scoreAngle*scoreDis + pre_Angle*dAngle1*scoreOrien*scoreDis)/(numTotal/2);
	                
	                der1[i+numTotal/2] += (pre_Orien*(-dOrien2)*scoreAngle*scoreDis + pre_Angle*dAngle2*scoreOrien*scoreDis + der_Dist*scoreOrien*scoreAngle)/(numTotal/2);
	                
	                der1[j+numTotal/2] += pre_Orien*dOrien2*scoreAngle*scoreDis/(numTotal/2);
	                
	                Vector ddfunc1, ddfunc2, ddfunc3, ddfunc4, ddfunc5;
	                ddfunc1 = pre_Orien*dOrien1*scoreAngle*scoreDis/(numTotal/2);
	                Tensor vv1(ddfunc1,dis1);
	                virial -= vv1;
	                
	                ddfunc2 = pre_Orien*dOrien2*scoreAngle*scoreDis/(numTotal/2);
	                Tensor vv2(ddfunc2,dis2);
	                virial -= vv2;
	                
	                ddfunc3 = pre_Angle*dAngle1*scoreOrien*scoreDis/(numTotal/2);
	                Tensor vv3(ddfunc3,dis1);
	                virial -= vv3;
	                
	                ddfunc4 = pre_Angle*dAngle2*scoreOrien*scoreDis/(numTotal/2);
	                Tensor vv4(ddfunc4,dist);
	                virial -= vv4;
	                
	                ddfunc5 = der_Dist*scoreOrien*scoreAngle/(numTotal/2);
	                Tensor vv5(ddfunc5,dist);
	                virial -= vv5;
                }
            }
        }
		for(int i = 0; i < der1.size(); i++)
        {
                derivatives.push_back(der1[i]);
        }

    }
    //Having calculated the cv, its derivative and the contribution to the virial you now transfer this information to the plumed core using the following three commands.
    for(unsigned int i=0;i<derivatives.size();i++)
    {
        setAtomsDerivatives(i,derivatives[i]);
    }
    //setBoxDerivatives(boxDerivatives);
    setBoxDerivatives(virial);
    //setBoxDerivativesNoPbc();
    setValue(cv_val);
}

}

}
