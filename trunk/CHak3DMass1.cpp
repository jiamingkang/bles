/*
 * CHak3DMass1.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak3DMass1.h"

CHak3DMass_1::CHak3DMass_1() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
    numNode=1;
	nodeSpace();
    volume = 1.0;
    mass = 0.0;
    connect=0;
    Trans=0;
	Ke=0;
	Me=0;
    massR=0;
    numR=0;
    T_loc=0;

}

CHak3DMass_1::~CHak3DMass_1() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
}

// compute transformation matrix and local coords
// Sandia Salinas RBE3 formulation
void CHak3DMass_1::makeTrans()
{
    int numConn = connect->getPos(); // number of nodes if the connection set
    Coord3D *length = new Coord3D [numConn];
    Coord3D ref; getNcrd(&ref); // mass node coordinate

    int i,j,ind,nd;
    for(i=0;i<numConn;i++)
    {
        nd = connect->nums[i]; // connected node number
        length[i] = gNodes[nd] - ref; // realtive coord to ref node
    }

    // make the S matrix
    int dofCon = 6*numConn;
    double *S = new double[6 * dofCon](); // 6n x 6
    double *W = new double[dofCon * dofCon](); // 6n x 6n
    for(i=0;i<numConn;i++)
    {
        ind = 36*i; // start of current block in S
        S[ind]=1.0; S[ind+7]=1.0; S[ind+14]=1.0;
        S[ind+21]=1.0; S[ind+28]=1.0; S[ind+35]=1.0; // ones on diagonal

        S[ind+4]=length[i].z; S[ind+5]=-length[i].y;
        S[ind+9]=-length[i].z; S[ind+11]=length[i].x;
        S[ind+15]=length[i].y; S[ind+16]=-length[i].x; // off diagonal coord terms

        // now for W matrix
        for(j=0;j<3;j++)
        {
            ind = i*6 + j; // dof number (translational dof)
            W[index2(ind,ind,dofCon)] = 1.0;
        }
    }

    // compute W S, 6n x 6
    double *WS = new double [6 * dofCon];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dofCon, 6, dofCon, 1.0, W, dofCon, S, 6, 0.0, WS, 6);

    // compute A = S^T W S
    double A[36]; // 6 x 6
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 6, 6, dofCon, 1.0, S, 6, WS, 6, 0.0, A, 6);

    // solve A^-1 (S^T W) to get the transfer matrix (A is symmetric as W is symmetric)
    j = dsy_slvLPK(6, dofCon, A, WS); // Note: WS sent as FORTRAN likes column major order

    // reudce the result (now in WS) to form the transfer matrix for translational dof
    int dofCon3 = 3*numConn;
    if(Trans){delete [] Trans;}
    Trans = new double [dofCon3*6]; // 6 rows, numConn x 3 columns

    // print out WS matrix
    /*for(i=0;i<dofCon;i++) // row
    {
        for(j=0;j<6;j++) // col
        {
            std::cout << "\nWS[" << i << "," << j << "] = " << WS[index2(i,j,6)];
        }
    }*/

    int ii,jj,rn,cn;
    for(i=0;i<numConn;i++)
    {
        j=6*i; // pointer to WS (row number)
        ind=3*i; // pointer to Trans (col number)

        // copy 9 entries from WS to Trans
        for(ii=0;ii<3;ii++)
        {
            rn = ii + j;
            cn = ii + ind;
            for(jj=0;jj<6;jj++)
            {
                Trans[index2(jj,cn,dofCon3)] = WS[index2(rn,jj,6)];
                //std::cout << "\nTrans[" << jj << "," << cn << "] = " << Trans[index2(jj,cn,dofCon3)];
            }
        }
    }

    // delete memory
    delete [] WS;
    delete [] S;
    delete [] W;
    delete [] length;
}

// compute transformation matrix and local coords
// Abaqus distributed coupling formulation
void CHak3DMass_1::makeTrans2()
{
    int i,j,ind,nd;
    int numConn = connect->getPos(); // number of nodes if the connection set
    Coord3D ref; getNcrd(&ref); // mass node coordinate
    double wgt = 1.0/(double)(numConn); // uniform weighting
    Coord3D x_bar; x_bar.x=0.0; x_bar.y=0.0; x_bar.z=0.0;

    // compute the weighted center
    for(i=0;i<numConn;i++)
    {
        nd = connect->nums[i]; // connected node number
        x_bar = x_bar + gNodes[nd]; // sum coords
    }
    x_bar = x_bar * wgt; // weighted center

    //compute realtive coords from nodes to x_bar
    Coord3D *rx = new Coord3D [numConn]; Coord3D rm;

    for(i=0;i<numConn;i++)
    {
        nd = connect->nums[i]; // connected node number
        rx[i] = gNodes[nd] - x_bar;
    }
    rm = ref - x_bar; // also for mass node

    // compute the intertia tensor (cor connection nodes)
    double *T = new double[9]();
    double r2;
    for(i=0;i<numConn;i++)
    {
        r2 = rx[i].x*rx[i].x + rx[i].y*rx[i].y + rx[i].z*rx[i].z; // dot product
        T[0] += r2; T[4] += r2; T[8] += r2; // diagonal

        // minus outer product (symmetric)
        T[0] -= rx[i].x*rx[i].x; T[1] -= rx[i].x*rx[i].y; T[2] -= rx[i].x*rx[i].z;
        T[4] -= rx[i].y*rx[i].y; T[5] -= rx[i].y*rx[i].z; T[8] -= rx[i].z*rx[i].z;
    }
    T[3] = T[1]; T[6] = T[2]; T[7] = T[5]; // symmetry

    // compute the R matrix (in column major order)
    int dofCon = 3*numConn;
    double *R = new double[3*dofCon](); // dofCon x 3

    for(i=0;i<numConn;i++)
    {
        ind=9*i; // start of current block
        R[ind+1] = rx[i].z; R[ind+2] = -rx[i].y;
        R[ind+3] = -rx[i].z; R[ind+5] = rx[i].x;
        R[ind+6] = rx[i].y; R[ind+7] = -rx[i].x;
    }

    // compute T^-1 R
    ind = dsy_slvLPK(3, dofCon, T, R); // Note: R in column major order

    // computr Rm matrix
    double Rm[9];
    Rm[0] = 0.0; Rm[1] = rm.z; Rm[2] = -rm.y;
    Rm[3] = -rm.z; Rm[4] = 0.0; Rm[5] = rm.x;
    Rm[6] = rm.y; Rm[7] = -rm.x; Rm[8] = 0.0;

    // multiply Rm x T^-1 x R
    double *Ttemp = new double[3*dofCon](); // 3 x dofCon
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, dofCon, 3, 1.0, Rm, 3, R, 3, 0.0, Ttemp, dofCon);

    // process the result (adding weight factor and additional values)
    if(Trans){delete [] Trans;}
    Trans = new double [dofCon*3]; // 3 rows, numConn x 3 columns
    for(i=0;i<3;i++) // row
    {
        nd=0;
        for(j=0;j<dofCon;j++) // col
        {
            ind = index2(i,j,dofCon); // indicator
            Trans[ind] = Ttemp[ind];
            if(i==nd){Trans[ind] += 1.0;}
            Trans[ind] *= wgt; //uniform weight
            nd++; if(nd==3){nd=0;} // indicate diagonal as we move accross the columns
            //std::cout << "\nTrans[" << i << "," << j << "] = " << Trans[ind];
        }
    }

    delete [] rx;
    delete [] T;
    delete [] R;
    delete [] Ttemp;
}

// function to add rigid mass connections
void CHak3DMass_1::setRigid(int nR, double *mR, int *ndR)
{
    int i,n;

    // re-set
    if(massR){delete [] massR;}
    if(nodeR){delete [] nodeR;}
    massR=0;
    nodeR=0;

    numR = nR;

    if(numR>0)
    {
        // copy data
        massR = new double [numR];
        nodeR = new int [numR];

        for(i=0;i<numR;i++)
        {
            massR[i] = mR[i];
            nodeR[i] = ndR[i];
            //std::cout << "\nmass " << massR[i] << " node num " << nodeR[i];
        }

        // compute the trasnfer matrices
        if(T_loc){delete [] T_loc;}
        T_loc = new double* [numR];

        Coord3D ref; getNcrd(&ref); // mass node coordinate
        Coord3D link;
        for(n=0;n<numR;n++)
        {
            T_loc[n] = new double [18](); // memory for 3 x 6 matrix
            for(i=0;i<3;i++){ T_loc[n][7*i]=1.0; } // 1's on diagonal, for left half

            // raotaion - translation in right half
            link = gNodes[nodeR[n]];
            T_loc[n][4] = link.z - ref.z; T_loc[n][9] = -T_loc[n][4];
            T_loc[n][15] = link.y - ref.y; T_loc[n][5] = -T_loc[n][15];
            T_loc[n][11] = link.x - ref.x; T_loc[n][16] = -T_loc[n][11];
        }
    }
}

// compute the mass matrix (mass transered matrix)
void CHak3DMass_1::makeMe()
{
    int i,j,k,ind;

    if(!Trans){makeTrans();} // ensure transfer matrix has been made

    int Me_ord = connect->getPos() * 3; // size of the mass matrix
    if(Me){delete [] Me;}
    float size = (Me_ord+1.0)*(Me_ord*0.5); // triangular storage space
    Me = new double [(int)size]();

    // first compute the local mass matrix
    double Me_loc[36]; // local mass matrix (6x6) to account for other masses connected by rigid elemen
    for(i=0;i<36;i++){Me_loc[i]=0.0;}

    for(i=0;i<3;i++){Me_loc[i*7] = mass;} // start with mass of "master" node

    double aux[36];
    for(i=0;i<numR;i++)
    {
        // matrix multiplication
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 6, 6, 3, 1.0, T_loc[i], 6, T_loc[i], 6, 0.0, aux, 6);

        // add to total
        for(j=0;j<36;j++){Me_loc[j] += aux[j]*massR[i];}
    }

    // now compute Me = Trans^T x Me_loc x Trans
    double *Mt = new double [Me_ord*6];
    double *Me_temp = new double [Me_ord*Me_ord];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, Me_ord, 6, 1.0, Me_loc, 6, Trans, Me_ord, 0.0, Mt, Me_ord);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Me_ord, Me_ord, 6, 1.0, Trans, Me_ord, Mt, Me_ord, 0.0, Me_temp, Me_ord);

    // then reduce to upper triangle
    for(i=0;i<Me_ord;i++) // row
    {
        for(j=i;j<Me_ord;j++) // col
        {
            ind = tri_ind(Me_ord,i,j); // upper triangle indicator
            k = index2(i,j,Me_ord); // full matrix indicator
            Me[ind] = Me_temp[k];
        }
    }

    delete [] Mt;
    delete [] Me_temp;
}

// add self-weight loading to a global loading vector
void CHak3DMass_1::addWeight(Coord3D vec, double mag, int *dof_ind, double *gload)
{
    if(!Trans){makeTrans();} // ensure transfer matrix has been made

    // first compute total element force (in each component)
    double force = mass * mag; // force = mass x accel
    double eload[6];
    eload[0] = vec.x * force;
    eload[1] = vec.y * force;
    eload[2] = vec.z * force; // force x vec
    eload[3] = 0.0; eload[4] = 0.0; eload[5] = 0.0; // zero moment force

    // add forces from rigid connections to eload
    int i,nd;
    double f_loc[3], f_tran[6];
    for(i=0;i<numR;i++)
    {
        force = massR[i] * mag;
        f_loc[0] = vec.x * force;
        f_loc[1] = vec.y * force;
        f_loc[2] = vec.z * force; // force x vec at connected node

        cblas_dgemv(CblasRowMajor, CblasTrans, 3, 6, 1.0, T_loc[i], 6, f_loc, 1, 0.0, f_tran, 1);
        for(nd=0;nd<6;nd++){ eload[nd] += f_tran[nd]; }
    }

    // expand force vector by multiplying by the transfer matrix
    int ord = connect->getPos() * 3; // size of the force vector
    double *geload = new double [ord]; // global force vector
    cblas_dgemv(CblasRowMajor, CblasTrans, 6, ord, 1.0, Trans, ord, eload, 1, 0.0, geload, 1);

    // add force x vec for each node into gloabl load vector (gload)
    int end = connect->getPos();
    int ind=0;
    for(i=0;i<end;i++)
    {
        //nd = getGlobal(i); // global node number
        nd = connect->nums[i];
        if(!dof_ind){nd*=3;} // assume 3 dof / node
        else{nd = dof_ind[nd];} // 1st degree of freedom (x)
        gload[nd++] += geload[ind++];
        gload[nd++] += geload[ind++]; // 2nd degree of freedom (y)
        gload[nd] += geload[ind++]; // 3rd degree of freedom (z)
    }

    delete [] geload;
}
