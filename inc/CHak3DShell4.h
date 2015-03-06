/*
 * CHak3DShell4.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK3DSHELL4_H_
#define CHAK3DSHELL4_H_

#include "inc/CHak3DElement.h"

class CHak3DShell_4: public CHak3DElement {
public:
	CHak3DShell_4();
	virtual ~CHak3DShell_4();
	double t;
	    double *Kstress;
	    Shell3* sub_shells[4]; // pointers to 4 sub elems

public:

	    // set thickness value
	    void setVals(double thk)
	    {
	        t=thk;
	        int i;
	        for(i=0;i<4;i++) { sub_shells[i]->setVals(thk); }
	    }
	    // delete stress stiffness matrix
	    void delKs()
	    {
	        if(Kstress)
	        {
	            delete [] Kstress; Kstress=0;
	            int i;
	            for(i=0;i<4;i++) { sub_shells[i]->delKs(); }
	        }
	    }
	    // set pointer to global node coordinate array
	    void setGNode(Coord3D *gn)
	    {
	        gNodes=gn; dof_ind=0;
	        int i;
	        for(i=0;i<4;i++) { sub_shells[i]->setGNode(gn); }
	    }
	    // set pointer to global node coordinate array & dof map
	    void setGNode(Coord3D *gn, int *dind)
	    {
	        gNodes=gn; dof_ind=dind;
	        int i;
	        for(i=0;i<4;i++) { sub_shells[i]->setGNode(gn,dind); }
	    }
	    // set material
	    void setMat(mat *m)
	    {
	        Mat=m;
	        int i;
	        for(i=0;i<4;i++) { sub_shells[i]->setMat(m); }
	    }
	    // return stess stiffness matrix
	    double* getKs(){return Kstress;}
	    // redefined addNode function (to add nodes to all sub elems)
	    void addNode(int loc,int glob)
	    {
	        nd[loc] = glob;

	        switch (loc)
			{
	            case 0:
	                sub_shells[0]->addNode(0,glob);
	                sub_shells[1]->addNode(0,glob);
	                sub_shells[2]->addNode(0,glob);
	                break;
				case 1:
	                sub_shells[0]->addNode(1,glob);
	                sub_shells[2]->addNode(1,glob);
	                sub_shells[3]->addNode(0,glob);
	                break;
				case 2:
	                sub_shells[0]->addNode(2,glob);
	                sub_shells[1]->addNode(1,glob);
	                sub_shells[3]->addNode(1,glob);
	                break;
				case 3:
	                sub_shells[1]->addNode(2,glob);
	                sub_shells[2]->addNode(2,glob);
	                sub_shells[3]->addNode(2,glob);
	                break;
			}
	    }
	    // function to make stiffness matrix
	    void makeKe()
	    {
	        int i, nt[3];
	        // make sub matrices 1st
	        for(i=0;i<4;i++) { sub_shells[i]->makeKe(); }

	        if(Ke){delete [] Ke;}
	        Ke = new double [300]();

	        // assemble sub matrices into total elem matrix
	        // first triangle is 0,1,2
			nt[0] = 0; nt[1] = 1; nt[2] = 2;
			addMatrix(3, 4, sub_shells[0]->getKe(), Ke, 6, nt);

			// second triangle is 0,2,3
			nt[1] = 2; nt[2] = 3;
			addMatrix(3, 4, sub_shells[1]->getKe(), Ke, 6, nt);

	        // third triangle is 0,1,3
			nt[1] = 1;
			addMatrix(3, 4, sub_shells[2]->getKe(), Ke, 6, nt);

	        // forth triangle is 1,2,3
			nt[0] = 1; nt[1] = 2;
			addMatrix(3, 4, sub_shells[3]->getKe(), Ke, 6, nt);

	        // now divide by 2
	        for(i=0;i<300;i++){Ke[i] *= 0.5;}
	    }
	    // function to make mass matrix (for sub elems)
	    void makeMe()
	    {
	        int i, nt[3];
	        for(i=0;i<4;i++) { sub_shells[i]->makeMe(); }

	        if(Me){delete [] Me;}
	        Me = new double [300]();

	        // assemble sub matrices into total elem matrix
	        // first triangle is 0,1,2
			nt[0] = 0; nt[1] = 1; nt[2] = 2;
			addMatrix(3, 4, sub_shells[0]->getMe(), Me, 6, nt);

			// second triangle is 0,2,3
			nt[1] = 2; nt[2] = 3;
			addMatrix(3, 4, sub_shells[1]->getMe(), Me, 6, nt);

	        // third triangle is 0,1,3
			nt[1] = 1;
			addMatrix(3, 4, sub_shells[2]->getMe(), Me, 6, nt);

	        // forth triangle is 1,2,3
			nt[0] = 1; nt[1] = 2;
			addMatrix(3, 4, sub_shells[3]->getMe(), Me, 6, nt);

	        // now divide by 2
	        for(i=0;i<300;i++){Me[i] *= 0.5;}
	    }
	    // function to make stress stiffness matrix (for sub elems)
	    void makeKs(double *disp)
	    {
	        int i, nt[3];
	        for(i=0;i<4;i++) { sub_shells[i]->makeKs(disp); }

	        if(Kstress){delete [] Kstress;}
	        Kstress = new double [300]();

	        // assemble sub matrices into total elem matrix
	        // first triangle is 0,1,2
			nt[0] = 0; nt[1] = 1; nt[2] = 2;
			addMatrix(3, 4, sub_shells[0]->getKs(), Kstress, 6, nt);

			// second triangle is 0,2,3
			nt[1] = 2; nt[2] = 3;
			addMatrix(3, 4, sub_shells[1]->getKs(), Kstress, 6, nt);

	        // third triangle is 0,1,3
			nt[1] = 1;
			addMatrix(3, 4, sub_shells[2]->getKs(), Kstress, 6, nt);

	        // forth triangle is 1,2,3
			nt[0] = 1; nt[1] = 2;
			addMatrix(3, 4, sub_shells[3]->getKs(), Kstress, 6, nt);

	        // now divide by 2
	        for(i=0;i<300;i++){Kstress[i] *= 0.5;}
	    }
	    // function to return von Mises stress (average of the four sub elems)
	    double vonMises(Coord3D p, double *disp)
	    {
	        int i; double svm = 0.0;
	        for(i=0;i<4;i++) { svm += sub_shells[i]->vonMises(p, disp); }
	        return svm*0.25;
	    }
	    // function to make volume (for sub elems & this elem)
	    void makeVol()
	    {
	        int i;
	        for(i=0;i<4;i++) { sub_shells[i]->makeVol(); }
	        volume = 0.5*( sub_shells[0]->getVol() + sub_shells[1]->getVol() );
	    }
	    // return thickness
	    double getThk(){return t;}
	    // compute Shell4 element thickness sensitivity
	    void elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
	                  double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf, int mode);
	    // compute Shell4 element thickness sensitivity w.r.t an eigenvalue
	    void EigSens(int numDual, int numEig, double **disp_prim, double **disp_dual,
	                 double *eig, double *wgt_case, double *elSens, int eNum, int mode);
	    // compute Shell4 flutter sensitivity - complex Eigenvalue
	    void CEigSens(int numDual, int numEig, dcompLPK **eig_flut, double **eig_vec,
	                  dcompLPK *eig, dcompLPK *wgt_case, double *elSens, int eNum);
	    // function to compute the contibution of an element to dks_du (buckling)
	    void DispSens(int numModes, double **disp_prim, double *fact, double *dk_du);
	    // function to compute the contibution of an element to dvm_du (stress)
	    void DispSensVM(int numModes, double **disp_prim, double *fact, double *dvm_du);
	    // add self-weight loading to a global loading vector
	    void addWeight(Coord3D vec, double mag, int *dof_ind, double *gload);
	    // function to update data due to a thickness change
	    void change_thk(double tin);
	    // output
		void KeOut() {
			int i;
			for(i=0;i<300;i++) {
				std::cout << "\n" << i << " - " << Ke[i];
			}
		}
		void MeOut() {
			int i;
			for(i=0;i<300;i++) {
				std::cout << "\n" << i << " - " << Me[i];
			}
		}


};

#endif /* CHAK3DSHELL4_H_ */
