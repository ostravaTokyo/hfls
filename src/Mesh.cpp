#include "Mesh.hpp"

using namespace std;


Mesh::Mesh()
{
    nElementsClst = 0;
    nPoints = 0;
    nSubClst = 0;
}

Mesh::~Mesh()
{

}



void Mesh::createMesh(){
    // - geometry setting ()
    double length[] = {4.0, 4.0, 4.0}; // corresponds to python benchmark
    //double length[] = {1.0, 1.0, 1.0};
    double radius = 1050.0;

    // - decomposition
    int nElSubXYZ[] = {5,5,5};
    int nSubXYZ[] = {4,4,4};


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                  GEOMETRY DEFINITION
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    double shift[3];
    int nElxyz_all[3];
    for (int i = 0; i < 3; i++){
        nElxyz_all[i] = nElSubXYZ[i] * nSubXYZ[i];
    }

    shift[0] = 0.5 * length[0];
    shift[1] = 0.5 * length[1];
    shift[2] =     0 * length[2];



    nSubClst = nSubXYZ[0] * nSubXYZ[1] * nSubXYZ[2];
    nElementsClst = nElxyz_all[0] * nElxyz_all[1] * nElxyz_all[2];
    nPoints = (nElxyz_all[0] + 1) * (nElxyz_all[1] + 1) * (nElxyz_all[2] + 1);


    elements.resize(nElementsClst);
    points.resize(nPoints);



    double dxyz[3];
    for (int i = 0 ; i < 3; i++){
        dxyz[i] = length[i] / nElxyz_all[i];
    }
    double _x,_y,_z;
    //  Points
    double weight;

    // corection of radius


    double Lz = length[2];
    double Lxy = 0.5 * sqrt(length[0] * length[0] + length[1] * length[1]);
    double radius_min = 0.5 *( Lz + (Lxy * Lxy) / Lz);
    bool flag_radius = true;
    if (radius < 1.1 * radius_min){
        radius = 1.1 *radius_min;
    }
    else if (radius > 100 * radius_min)
    {
        flag_radius = false;
    }





    int cnt = 0;
    for (int kk = 0; kk <  nElxyz_all[2] + 1; kk++){
        for (int jj = 0; jj <  nElxyz_all[1] + 1; jj++){
            for (int ii = 0; ii <  nElxyz_all[0] + 1; ii++){
                _x = ii * dxyz[0] - shift[0];
                _y = jj * dxyz[1] - shift[1];
                _z = kk * dxyz[2];
                if (flag_radius){
                    weight = (length[2] - _z) / length[2];
                    _z = (-sqrt(pow(radius,2) -
                            (pow(_x,2) + pow(_y,2))) + radius) * weight + _z;
                }
                _x += 0.5 * length[0];
                _y += 0.5 * length[1];
                _z += 0.0 * length[2];
                points[cnt].x =_x;
                points[cnt].y =_y;
                points[cnt].z =_z;
                cnt++;
            }
        }
    }

    int ttt[8], _first_set[4];
    vector < int > tmp_vec (ttt, ttt + sizeof(ttt)/sizeof(int));
    int nxy= (nElxyz_all[0] + 1) * (nElxyz_all[1] + 1);
    cnt = 0;
    for (int kk = 0; kk <  nElxyz_all[2]; kk++){
        for (int jj = 0; jj <  nElxyz_all[1]; jj++){
            for (int ii = 0; ii <  nElxyz_all[0]; ii++){
                _first_set[0] = ii + 0 + (jj + 0) * (nElxyz_all[0] + 1);
                _first_set[1] = ii + 1 + (jj + 0) * (nElxyz_all[0] + 1);
                _first_set[2] = ii + 1 + (jj + 1) * (nElxyz_all[0] + 1);
                _first_set[3] = ii + 0 + (jj + 1) * (nElxyz_all[0] + 1);
                for (int ll = 0; ll < 4; ll++){
                    tmp_vec[ll]     = _first_set[ll] + kk * nxy;
                    tmp_vec[ll+4]   = _first_set[ll] + (kk + 1) * nxy;
                }
                for (int ll = 0; ll < 8; ll ++){
                    elements[cnt].ind[ll] = tmp_vec[ll];
                }
                cnt++;
            }
        }
    }

    //

    //
    //---------------------------
    //
    cnt = 0;
    int currentIdOfMat = 0;
    for (int KK = 0; KK <  nSubXYZ[2]; KK++){
        for (int JJ = 0; JJ <  nSubXYZ[1]; JJ++){
            for (int II = 0; II <  nSubXYZ[0]; II++){
                for (int kk = KK * nElSubXYZ[2] ; kk <  (KK + 1) * nElSubXYZ[2]; kk++){
                    for (int jj = JJ * nElSubXYZ[1]; jj <  (JJ + 1) * nElSubXYZ[1]; jj++){
                        for (int ii = II * nElSubXYZ[0]; ii <  (II + 1) * nElSubXYZ[0]; ii++){
                            cnt = ii + jj * nElxyz_all[0] + kk * (nElxyz_all[0] * nElxyz_all[1]);
                            elements[cnt].PartitionId =  currentIdOfMat;;
                        }
                    }
                }
                currentIdOfMat++;
            }
        }
    }
/*
    */
}

//void Mesh::feti_ddm_cluster(){
//
//
//    if (nSubClst == 0){
//        //nSubClst = nSubXYZ[0] * nSubXYZ[1] + nSubXYZ[2];
//
//}
