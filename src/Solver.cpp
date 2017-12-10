// ab
#include "Solver.hpp"



Solver::Solver()
{

}

void Solver::pcpg(map <string,string> options2, Cluster &cluster){
      clock_t begin = clock();

    double eps_iter     = atof(options2["eps_iter"].c_str());
    double max_iter     = atoi(options2["max_iter"].c_str());

    double gPz, gPz_prev, wFw, rho, gamma, norm_gPz0;
    Vector g0, d_rhs, e, iGTG_e, lambda, z, Pz;
    Vector Fw, Pg, g, w, w_prev;
    Vector beta, alpha;

    vector < Vector > xx, yy;

    int nSubClst = cluster.get_nSubClst();
    xx.resize(nSubClst);
    yy.resize(nSubClst);

    xx[0].label = "test";
    cluster.mult_Kplus_f(cluster.rhs,xx);

    cluster.mult_Bf(xx,d_rhs);

    // e = Rt * f
    cluster.mult_RfT(cluster.rhs,e);

    // lambda0 = G * inv(GtG) * e
    iGTG_e.mat_mult_dense(cluster.invGfTGf,"N",e,"N");
    cluster.mult_Gf(iGTG_e, lambda);


    // F * lambda0
    cluster.mult_Ff(lambda,g0);
    // g0 = F * lambda0 - d_rhs
    g0.add(d_rhs,-1);

    g = g0;

    // Pg0
    cluster.Projection(g,Pg,alpha);
    cluster.Preconditioning(Pg,z);
    cluster.Projection(z,Pz,beta);
    gPz = Matrix::dot(g,Pz);

    double norm_g0Pg0 = sqrt(Matrix::dot(g,Pg));

    norm_gPz0 = sqrt(gPz);

    printf("\n|g0Pg0| = %3.9e  \n", norm_g0Pg0);
    printf(  "|gPz0|  = %3.9e  \n\n", norm_gPz0);
    w = Pz;

    printf("=======================\n");
    printf(" it.\t||gradient||\n");
    printf("=======================\n");

    for (int it = 0; it < max_iter; it++){

        printf("%4d\t%3.9e \n",it + 1, sqrt(gPz) / norm_gPz0);
        if (sqrt(gPz) < eps_iter * norm_gPz0)
            break;

        cluster.mult_Ff(w,Fw);
        wFw = Matrix::dot(w,Fw);
        rho = -gPz / wFw;

        lambda.add(w,rho);
        g.add(Fw,rho);

        cluster.Projection(g,Pg,alpha);
        cluster.Preconditioning(Pg,z);
        cluster.Projection(z,Pz,beta);

        gPz_prev = gPz;
        gPz = Matrix::dot(g,Pz);

        gamma = gPz / gPz_prev;

        w_prev = w;
        w = Pz;
        w.add(w_prev,gamma);

        if (options2["vtkWithinIter"].compare("true") == 0)
            cluster.printVTK(yy, xx, lambda, alpha, it);
    }
    clock_t end = clock();
    time_solver = double(end - begin) / CLOCKS_PER_SEC;


    // final solution (last parameter -1 avoids numbering of vtk file)
    printf("Solver time:  %3.1f s.\n",time_solver);
    printf("Total time:   %3.1f s.\n",time_total);

    cluster.printVTK(yy, xx, lambda, alpha, -1);

}
