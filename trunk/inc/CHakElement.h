/*
 * CHakElement.h
 *
 *  Created on: 22 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAKELEMENT_H_
#define CHAKELEMENT_H_

class CHakElement {
public:
	CHakElement();
	virtual ~CHakElement();
protected:
	int numNode;
	int *nd;	// array of node numbers
	mat *Mat;	// pointer to element material class
	double *Ke;	// array for stiffness matrix
	double *Me;	// array for mass martix
    double volume; // volume of the element (useful for self-weight loading)
	Coord3D *gNodes; // pointer to the gloabl node array
    int *dof_ind; // pointer to a dof map for the mesh that the elem belongs to
public:
	virtual void addNode(int loc,int glob){nd[loc] = glob;} // add node number to element
	int getGlobal(int loc) {return(nd[loc]);}		// get global node num from local num
	int* getNodes(){return(nd);}					// get pointer to element nodes
    int getNum(){return(numNode);}
	void nodeSpace(){nd=new int [numNode];}	// set storage for global node num array
	virtual void setGNode(Coord3D *gn){gNodes=gn; dof_ind=0;}	// set pointer to global node coordinate array
    virtual void setGNode(Coord3D *gn, int *dind){gNodes=gn; dof_ind=dind;}	// set pointer to global node coordinate array & dof map
	virtual void setMat(mat *m){Mat=m;}				// set material
    mat* getMat(){return Mat;}
	// virtual destructor
	virtual ~elem(){if(nd){delete []nd;}}
	// read in node coordinates
	void getNcrd(Coord3D *np) {
		int i,j;
		for(i=0;i<numNode;i++) {
			j = getGlobal(i);
			np[i].x = gNodes[j].x;
			np[i].y = gNodes[j].y;
			np[i].z = gNodes[j].z;
		} }
	// create stiffness matrix
	virtual void makeKe()=0;
	// create mass matrix
	virtual void makeMe()=0;
    // compute volume
    virtual void makeVol()=0;
    // function to clear matrices
    void delMat(){
        if(Ke){delete [] Ke;Ke=0;}
		if(Me){delete [] Me;Me=0;}}
	// get matrix pointers
	double* getKe(){return(Ke);}
	double* getMe(){return(Me);}
    double getVol(){return volume;}
	virtual void KeOut()=0;
	virtual void MeOut()=0;

};

#endif /* CHAKELEMENT_H_ */
