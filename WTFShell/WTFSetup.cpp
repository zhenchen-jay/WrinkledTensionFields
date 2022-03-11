#include "WTFSetup.h"
#include "../GeometryDerivatives.h"
#include "../MeshConnectivity.h"

void WTF::WTFSetup::buildRestFundamentalForms()
{
    int nfaces = restF.rows();
    abars.resize(nfaces);
    bbars.resize(nfaces);

    MeshConnectivity restMesh(restF);
    sff->initializeExtraDOFs(restEdgeDOFs, restMesh, restV);

    for (int i = 0; i < nfaces; i++)
    {
        abars[i] = firstFundamentalForm(restMesh, restV, i, NULL, NULL);
        if (restFlat)
            bbars[i].setZero();
        else
            bbars[i] = sff->secondFundamentalForm(restMesh, restV, restEdgeDOFs, i, NULL, NULL);
    }
}

void WTF::WTFSetup::buildQuadraturePoints()
{
    if (quadNum == 1)
    {
        quadPoints.clear();
        QuadraturePoints point;

        point.u = 1.0 / 3;
        point.v = 1.0 / 3;
        point.weight = 1.0;
        quadPoints.push_back(point);
    }
    else if (quadNum == 3)
    {
        quadPoints.clear();
        QuadraturePoints point;

        double x, pos[3];
        x = 0.16666666666666666667;
        point.weight = 0.333333333333333333;
        pos[0] = x;
        pos[1] = x;
        pos[2] = 1 - 2 * x;
        for (int i = 0; i < 3; i++)
        {
            point.u = pos[i];
            point.v = pos[(i + 1) % 3];
            quadPoints.push_back(point);
        }
    }
    else if (quadNum == 6)
    {
        quadPoints.clear();
        QuadraturePoints point;

        double x, pos[3];
        x = 0.091576213509771;
        point.weight = 0.109951743655322;
        pos[0] = x;
        pos[1] = x;
        pos[2] = 1 - 2 * x;
        for (int i = 0; i < 3; i++)
        {
            point.u = pos[i];
            point.v = pos[(i + 1) % 3];
            quadPoints.push_back(point);
        }

        x = 0.445948490915965;
        point.weight = 0.223381589678011;
        pos[0] = x;
        pos[1] = x;
        pos[2] = 1 - 2 * x;
        for (int i = 0; i < 3; i++)
        {
            point.u = pos[i];
            point.v = pos[(i + 1) % 3];
            quadPoints.push_back(point);
        }
    }
    else if (quadNum == 16)
    {
        quadPoints.clear();
        QuadraturePoints point;

        point.u = 1.0 / 3;
        point.v = 1.0 / 3;
        point.weight = 0.1443156076777871682510911104890646;
        quadPoints.push_back(point);

        double x, pos[3];
        x = 0.1705693077517602066222935014914645;
        point.weight = 0.1032173705347182502817915502921290;
        pos[0] = x;
        pos[1] = x;
        pos[2] = 1 - 2 * x;
        for (int i = 0; i < 3; i++)
        {
            point.u = pos[i];
            point.v = pos[(i + 1) % 3];
            quadPoints.push_back(point);
        }

        x = 0.0505472283170309754584235505965989;
        point.weight = 0.0324584976231980803109259283417806;
        pos[0] = x;
        pos[1] = x;
        pos[2] = 1 - 2 * x;
        for (int i = 0; i < 3; i++)
        {
            point.u = pos[i];
            point.v = pos[(i + 1) % 3];
            quadPoints.push_back(point);
        }

        x = 0.459292588292723156028815514494169;
        point.weight = 0.0950916342672846247938961043885843;
        pos[0] = x;
        pos[1] = x;
        pos[2] = 1 - 2 * x;
        for (int i = 0; i < 3; i++)
        {
            point.u = pos[i];
            point.v = pos[(i + 1) % 3];
            quadPoints.push_back(point);
        }

        x = 0.263112829634638113421785786284643;
        double y = 0.008394777409957605337213834539294;
        point.weight = 0.0272303141744349942648446900739089;
        pos[0] = x;
        pos[1] = y;
        pos[2] = 1 - x - y;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 1; j < 3; j++)
            {
                point.u = pos[i];
                point.v = pos[(i + j) % 3];
                quadPoints.push_back(point);
            }
        }
    }
    else
    {
        std::cout << "Please specify the order of quadrature points: 1, 3, 6, 16. Set default to 3" << std::endl;
        quadPoints.clear();
        QuadraturePoints point;

        double x, pos[3];
        x = 0.16666666666666666667;
        point.weight = 0.333333333333333333;
        pos[0] = x;
        pos[1] = x;
        pos[2] = 1 - 2 * x;
        for (int i = 0; i < 3; i++)
        {
            point.u = pos[i];
            point.v = pos[(i + 1) % 3];
            quadPoints.push_back(point);
        }

    }

}

