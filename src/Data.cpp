#include "Data.hpp"

Data::Data()
{
}


Data::~Data()
{

}


void Data::fe_assemb_local_K_f(Mesh &mesh){

    Point coord[8];
    local_K_f_clust.resize(mesh.nElementsClst);

    double E = mesh.material.young_modulus;
    double mu = mesh.material.poissons_ratio;

    for (int i = 0; i < mesh.nElementsClst; i++){
        Element &i_elem =  mesh.elements[i];
        for (int j = 0; j < 8; j++){
            int iii = i_elem.ind[j];
            //
            coord[j].x = mesh.points[iii].x;
            coord[j].y = mesh.points[iii].y;
            coord[j].z = mesh.points[iii].z;
            //
            local_K_f_clust[i].ieq[j + 0] = 3 * iii + 0;
            local_K_f_clust[i].ieq[j + 8] = 3 * iii + 1;
            local_K_f_clust[i].ieq[j + 16] = 3 * iii + 2;
        }
        stf_mtrx_solid45(local_K_f_clust[i], coord, E, mu);
    }

//    delete [] coord;

}


int compareInt(const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

void Data::create_analytic_ker_K(Mesh &mesh, vector <Matrix> &R_){

    int nSubClst = mesh.nSubClst;


    for (int d = 0; d < nSubClst; d++){
        Matrix &R = R_[d];
        int n = l2g[d].size();
        R.zero_dense(n,6);
        double x,y,z;
        int n_nods = int(n / 3);
//        cout<< "n       =    "<< n << endl;
//        cout<< "n_nods  =    "<< n_nods << endl;
        for (int i = 0 ; i < n_nods; i++){

            int in_glob = int ( l2g[d][3*i] / 3);

            x = mesh.points[in_glob].x;
            y = mesh.points[in_glob].y;
            z = mesh.points[in_glob].z;


            R.dense[3 * i + 0 + n * 0] = 1;
            R.dense[3 * i + 1 + n * 1] = 1;
            R.dense[3 * i + 2 + n * 2] = 1;

            R.dense[3 * i + 0 + n * 3] = -y;
            R.dense[3 * i + 1 + n * 3] =  x;

            R.dense[3 * i + 0 + n * 4] = -z;
            R.dense[3 * i + 2 + n * 4] =  x;

            R.dense[3 * i + 1 + n * 5] = -z;
            R.dense[3 * i + 2 + n * 5] =  y;
        }

        for (int j = 0; j < R.n_col; j++){
            double norm_2 = 0;
            for (int i = 0; i < R.n_row; i++){
               norm_2 += R.dense[i + n * j] *  R.dense[i + n * j];
            }
            norm_2 = sqrt(norm_2);
            for (int i = 0; i < R.n_row; i++){
               R.dense[i + n * j] /= norm_2;
            }
        }

    }
}


void Data::feti_symbolic(Mesh &mesh, vector <Matrix> &K_)
{

    //
    selectorOfElemPartitId.resize(mesh.nSubClst);
    std::vector<int> l2g_dim(mesh.nSubClst, 0);

    for (int i = 0 ; i < mesh.nElementsClst ; i++){
        l2g_dim[mesh.elements[i].PartitionId] +=  local_K_f_clust[i].nDOF;
        selectorOfElemPartitId[mesh.elements[i].PartitionId].push_back(i);
    }

   //
    l2g.resize(mesh.nSubClst);
    g2l.resize(mesh.nSubClst);
    //
    for (int d = 0 ; d < mesh.nSubClst ; d++){
        l2g[d].resize(l2g_dim[d]);
        int cnt = 0;
        for (int i = 0 ; i < selectorOfElemPartitId[d].size() ; i++){
            local_K_f &i_loc_K_f = local_K_f_clust[selectorOfElemPartitId[d][i]];
            for (int k = 0 ; k < i_loc_K_f.nDOF; k++){
                l2g[d][k + cnt] = i_loc_K_f.ieq[k];
            }
            cnt += i_loc_K_f.nDOF;
        }
        qsort(&(l2g[d][0]), l2g[d].size(), sizeof (int), compareInt);
        std::vector<int>::iterator itv;
        itv = std::unique (l2g[d].begin(), l2g[d].end());
        l2g[d].resize( std::distance(l2g[d].begin(),itv) );

        for (unsigned int i = 0; i < l2g[d].size(); i++){
            g2l[d].insert ( pair<int,int>(l2g[d][i],i) );
        }


        /* renumbering stiff mat from global to local (subdom.) numbering */
        for (int i = 0 ; i < selectorOfElemPartitId[d].size() ; i++){
            for (int k = 0 ; k < local_K_f_clust[selectorOfElemPartitId[d][i]].nDOF; k++){
                local_K_f &i_loc_K_f = local_K_f_clust[selectorOfElemPartitId[d][i]];
                i_loc_K_f.ieq[k] = g2l[d][i_loc_K_f.ieq[k]];
            }
        }

    }



    // variable 'bool symMatrix' prepared for nonsymm. cases

    vector<int>::iterator itv;
    for (int d = 0 ; d < mesh.nSubClst; d++){
        bool symMatrix=true;
        int *indK;
        vector < vector < int > > forCSRformat;
        //forCSRformat.resize(fem->domain->neqSub,vector<int>(0));
        forCSRformat.resize(l2g[d].size(),vector <int> (0));
    //
        for (int i = 0 ; i < selectorOfElemPartitId[d].size() ; i++){
            local_K_f &i_loc_K_f = local_K_f_clust[selectorOfElemPartitId[d][i]];
            int elemDOFs = i_loc_K_f.nDOF;
            indK = i_loc_K_f.ieq;
            for (int j=0;j<elemDOFs;j++){
              for (int k=0;k<elemDOFs;k++){
                if (symMatrix && (indK[k]>=indK[j]) || !symMatrix){
                  forCSRformat[indK[j]].push_back(indK[k]);
                }
              }
            }
        }
        int nnz_K = 0;
        for (int i = 0;i<forCSRformat.size();i++){
          sort(forCSRformat[i].begin(),forCSRformat[i].end());
          itv = std::unique (forCSRformat[i].begin(), forCSRformat[i].end());
          forCSRformat[i].resize( std::distance(forCSRformat[i].begin(),itv) );
          nnz_K+=forCSRformat[i].size();
        }



        int neqSub = l2g[d].size();
        Matrix &_K = K_[d];
        _K.label = "stiffness matrix";
        _K.order_number = d;

        _K.j_col.resize(nnz_K);
        _K.val.resize(nnz_K);
        for (int i = 0; i < nnz_K; i++)
            _K.val[i] = 0;
        _K.i_ptr.resize(neqSub + 1);
        int cnt_it  = 0;
        _K.i_ptr[0]  = 0;
        for (int i=0;i<forCSRformat.size();i++){
          for (vector<int>::iterator it1 = forCSRformat[i].begin();it1!=forCSRformat[i].end();it1++){
            _K.j_col[cnt_it]	= *it1;
            cnt_it++;
            _K.i_ptr[i+1]=cnt_it;
          }
        }
        _K.format = 1;
        _K.symmetric = 2;
        _K.nnz = nnz_K;
        _K.n_row_cmprs = neqSub;
        _K.n_row = neqSub;
        _K.n_col = neqSub;
        _K.l2g_i_coo.resize(neqSub);
        for (int i = 0 ; i < neqSub; i++){
            _K.l2g_i_coo[i] = i;
        }

    }

}


void Data::feti_numeric_element(Matrix &Ksub, Vector & rhs_sub, local_K_f &Kelem){
      //

    double K_ij;
    int i_ind, j_ind;
    int nDOF = Kelem.nDOF;
    for (int i = 0; i < nDOF; i++) {

        rhs_sub.dense[Kelem.ieq[i]] += Kelem.val_f[i];

        for (int j = i; j < nDOF; j++) {
            j_ind = Kelem.ieq[j];
            K_ij = Kelem.val_K[i * nDOF + j];
            if (j_ind < Kelem.ieq[i]) {
                i_ind = Kelem.ieq[j];
                j_ind = Kelem.ieq[i];
            }
            else  {
                i_ind = Kelem.ieq[i];
            }
            for (int k = Ksub.i_ptr[i_ind]; k < Ksub.i_ptr[i_ind + 1]; k++) {
                if ((Ksub.j_col[k]) == j_ind)  {
                    Ksub.val[k] += K_ij;
                break;
                } // if
            } // loop: k
        } // loop: j
    } // loop: i
}



void Data::feti_numeric(Mesh &mesh, vector <Matrix> & K_,  vector <Vector>& rhs_)
{
    for (int d = 0 ; d < mesh.nSubClst; d++){
      //        rhs_[d].zero_dense(K_[d].n_row_cmprs,1);
      rhs_[d].zero_dense(K_[d].n_row_cmprs);
        for (int i = 0 ; i < selectorOfElemPartitId[d].size() ; i++){
            local_K_f &i_loc_K_f = local_K_f_clust[selectorOfElemPartitId[d][i]];
            feti_numeric_element(K_[d],rhs_[d], i_loc_K_f);
        }
    }
}


double Data::inverse_matrix_3x3(double *A, double *iA) {

    double determinant = +A[0] * (A[4] * A[8] - A[5] * A[7])
            - A[3] * (A[1] * A[8] - A[7] * A[2])
            + A[6] * (A[1] * A[5] - A[4] * A[2]);
    double invdet = 1. / determinant;

    iA[0] = (A[4] * A[8] - A[5] * A[7]) * invdet;
    iA[1] = -(A[3] * A[8] - A[6] * A[5]) * invdet;
    iA[2] = (A[3] * A[7] - A[6] * A[4]) * invdet;
    iA[3] = -(A[1] * A[8] - A[7] * A[2]) * invdet;
    iA[4] = (A[0] * A[8] - A[6] * A[2]) * invdet;
    iA[5] = -(A[0] * A[7] - A[1] * A[6]) * invdet;
    iA[6] = (A[1] * A[5] - A[2] * A[4]) * invdet;
    iA[7] = -(A[0] * A[5] - A[2] * A[3]) * invdet;
    iA[8] = (A[0] * A[4] - A[1] * A[3]) * invdet;

    return determinant;
}


void Data::stf_mtrx_solid45(local_K_f & i_local_K_f, Point *coordinate, double E, double mu){




//    for (int i = 0; i < 8; i++){
//       cout << coordinate[i].x << " " ;
//       cout << coordinate[i].y << " " ;
//       cout << coordinate[i].z << "\n";
//    }

    /* hexahedron
     * local element matrix
     * hexahedron (in ANSYS - solid 45)
     * 8 nodes,3 degrees of freedom per node
     * 24x24, 576 entries
     */
    i_local_K_f.nDOF = 24;
    double density = 2.;

    //double acceleration[] ={0.0,0.0, 98.10};
    double acceleration[] ={0.0, 0 * 0.09810,98.10};

    double r, s, t;
    double Gt[6 * 24];
    double N[8], dNr[8], dNs[8], dNt[8]; // dNx[8], dNy[8], dNz[8]
    double Jmat[9], iJmat[9];
    double CG_tmp[6 * 24];
    double dNx_j, dNy_j, dNz_j;

    memset(i_local_K_f.val_f,0,24*sizeof(double));
    memset(i_local_K_f.val_K,0,576*sizeof(double));
    memset(Gt,0,144*sizeof(double));
    memset(N,0,8*sizeof(double));


    double C [36];
    double gaussPoints[24];
//    double E = 10000;
    //double E = 2.1e9;  // if this mat. is set up,
    //                      Dissection has problem with accuracy
    //                      during the factorization of Ac = [ Fc  Gc ]
    //                                                       [ Gct  O ] 
 //   double mu = 0.3;
    double const1 = E / ((1. + mu) * (1. - 2. * mu));
    double mu2 = 1. - mu;
    double mu3 = 0.5 - mu;
    double pt = 0.577350269189626;
    //
    memset(&(gaussPoints[0]),0, 24 * sizeof(double));
    memset(&(C[0]),0, 36* sizeof(double));
    //
    C[0] = mu2 * const1;
    C[6] = mu * const1;
    C[12] = mu * const1;
    C[1] = mu * const1;
    C[7] = mu2 * const1;
    C[13] = mu * const1;
    C[2] = mu * const1;
    C[8] = mu * const1;
    C[14] = mu2 * const1;
    C[21] = mu3 * const1;
    C[28] = mu3 * const1;
    C[21] = mu3 * const1;
    C[35] = mu3 * const1;
    //
    gaussPoints[0] = -pt; // 0
    gaussPoints[1] = -pt;
    gaussPoints[2] = -pt;
    gaussPoints[3] = -pt;
    gaussPoints[4] = pt;
    gaussPoints[5] = pt;
    gaussPoints[6] = pt;
    gaussPoints[7] = pt;
    gaussPoints[8] = -pt; // 1
    gaussPoints[9] = -pt;
    gaussPoints[10] = pt;
    gaussPoints[11] = pt;
    gaussPoints[12] = -pt;
    gaussPoints[13] = -pt;
    gaussPoints[14] = pt;
    gaussPoints[15] = pt;
    gaussPoints[16] = -pt; // 2
    gaussPoints[17] = pt;
    gaussPoints[18] = -pt;
    gaussPoints[19] = pt;
    gaussPoints[20] = -pt;
    gaussPoints[21] = pt;
    gaussPoints[22] = -pt;
    gaussPoints[23] = pt;



    for (int i = 0; i < 8; i++) {
        //
        r = gaussPoints[i];
        s = gaussPoints[i + 8];
        t = gaussPoints[i + 16];
        //
        dNr[0] = 0.125 * -(1. - s) * (1. - t);
        dNr[1] = 0.125 * (1. - s) * (1. - t);
        dNr[2] = 0.125 * (1. + s) * (1. - t);
        dNr[3] = 0.125 * -(1. + s) * (1. - t);
        dNr[4] = 0.125 * -(1. - s) * (1. + t);
        dNr[5] = 0.125 * (1. - s) * (1. + t);
        dNr[6] = 0.125 * (1. + s) * (1. + t);
        dNr[7] = 0.125 * -(1. + s) * (1. + t);
        //
        dNs[0] = 0.125 * -(1. - r) * (1. - t);
        dNs[1] = 0.125 * -(1. + r) * (1. - t);
        dNs[2] = 0.125 * (1. + r) * (1. - t);
        dNs[3] = 0.125 * (1. - r) * (1. - t);
        dNs[4] = 0.125 * -(1. - r) * (1. + t);
        dNs[5] = 0.125 * -(1. + r) * (1. + t);
        dNs[6] = 0.125 * (1. + r) * (1. + t);
        dNs[7] = 0.125 * (1. - r) * (1. + t);
        //
        dNt[0] = 0.125 * -(1. - r) * (1. - s);
        dNt[1] = 0.125 * -(1. + r) * (1. - s);
        dNt[2] = 0.125 * -(1. + r) * (1. + s);
        dNt[3] = 0.125 * -(1. - r) * (1. + s);
        dNt[4] = 0.125 * (1. - r) * (1. - s);
        dNt[5] = 0.125 * (1. + r) * (1. - s);
        dNt[6] = 0.125 * (1. + r) * (1. + s);
        dNt[7] = 0.125 * (1. - r) * (1. + s);
        //
        N[0] = 0.125 * (1. - r) * (1. - s) * (1. - t);
        N[1] = 0.125 * (r + 1.) * (1. - s) * (1. - t);
        N[2] = 0.125 * (r + 1.) * (s + 1.) * (1. - t);
        N[3] = 0.125 * (1. - r) * (s + 1.) * (1. - t);
        N[4] = 0.125 * (1. - r) * (1. - s) * (t + 1.);
        N[5] = 0.125 * (r + 1.) * (1. - s) * (t + 1.);
        N[6] = 0.125 * (r + 1.) * (s + 1.) * (t + 1.);
        N[7] = 0.125 * (1. - r) * (s + 1.) * (t + 1.);
        //
        for (int j = 0; j < 9; j++) {
            Jmat[j] = 0.;
            iJmat[j] = 0.;
        }
        //
        for (int j = 0; j < 8; j++) {
            Jmat[0] += (dNr[j] * coordinate[j].x);
            Jmat[1] += (dNs[j] * coordinate[j].x);
            Jmat[2] += (dNt[j] * coordinate[j].x);
            Jmat[3] += (dNr[j] * coordinate[j].y);
            Jmat[4] += (dNs[j] * coordinate[j].y);
            Jmat[5] += (dNt[j] * coordinate[j].y);
            Jmat[6] += (dNr[j] * coordinate[j].z);
            Jmat[7] += (dNs[j] * coordinate[j].z);
            Jmat[8] += (dNt[j] * coordinate[j].z);
        }
        //
        double detJ = inverse_matrix_3x3(Jmat, iJmat);
        //
        for (int j = 0; j < 8; j++) {
            dNx_j = iJmat[0] * dNr[j] + iJmat[3] * dNs[j] + iJmat[6] * dNt[j];
            dNy_j = iJmat[1] * dNr[j] + iJmat[4] * dNs[j] + iJmat[7] * dNt[j];
            dNz_j = iJmat[2] * dNr[j] + iJmat[5] * dNs[j] + iJmat[8] * dNt[j];
            Gt[j] = dNx_j; // x
            Gt[j + 32] = dNy_j; // y
            Gt[j + 64] = dNz_j; // z
            Gt[j + 72] = dNy_j; // y-x
            Gt[j + 80] = dNx_j; // y-x
            Gt[j + 104] = dNz_j; // y-z
            Gt[j + 112] = dNy_j; // y-z
            Gt[j + 120] = dNz_j; // x-z
            Gt[j + 136] = dNx_j; // x-z
        }
        //
        memset(CG_tmp,0,144*sizeof(double));
        //
        for (int I = 0; I < 6; I++) {
          for (int J = 0; J < 24; J++) {
            CG_tmp[I * 24 + J] = 0.;
            for (int K = 0; K < 6; K++) {
              CG_tmp[I * 24 + J] += (C[I + K * 6] * Gt[24 * K + J]);
            }
          }
        }
        for (int I = 0; I < 24; I++) {
          for (int J = 0; J < 24; J++) {
            for (int K = 0; K < 6; K++) {
              i_local_K_f.val_K[I + J * 24] += (Gt[I + K * 24]
                  * CG_tmp[K * 24 + J] * detJ);
            }
          }
        }
        //
        for (int j = 0; j < 8; j++) {
          i_local_K_f.val_f[j] += acceleration[0] * density * N[j] * detJ;
          i_local_K_f.val_f[j + 8] += acceleration[1] * density * N[j] * detJ;
          i_local_K_f.val_f[j + 16] += acceleration[2] * density * N[j] * detJ;
        }
    }
}


