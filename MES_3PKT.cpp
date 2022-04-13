#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double kodt = 25.; //przewodnosc
double DOKLADNOSC = 1.e-10;
double H = 0.1, B = 0.1;
const int nH = 4, nB = 4;
double tempOt = 1200.; //temperatura otoczenia
double alpha = 300.;
double c = 700.; // cieplo wlasciwe 
double ro = 7800.; // gestosc
double t0 = 100.; //temperatura poczatkowa
double t = 500.; //czas procesu
double dt = 50.; //krok czasowy

struct Pkt {
    double pkt[4][3] = {};
    double pktCal[9][4] = {
        //x wagax y wagay
        {-sqrt(3. / 5.), 5. / 9., -sqrt(3. / 5.), 5. / 9.},
        {0., 8. / 9., -sqrt(3. / 5.), 5. / 9.},
        {sqrt(3. / 5.), 5. / 9., -sqrt(3. / 5.), 5. / 9.},

        {-sqrt(3. / 5.), 5. / 9., 0., 8. / 9.},
        {0., 8. / 9., 0., 8. / 9.},
        {sqrt(3. / 5.), 5. / 9., 0., 8. / 9.},

        {-sqrt(3. / 5.), 5. / 9., sqrt(3. / 5.), 5. / 9.},
        {0., 8. / 9., sqrt(3. / 5.), 5. / 9.},
        {sqrt(3. / 5.), 5. / 9., sqrt(3. / 5.), 5. / 9.},
    };
    double funkcjeKsztaltu[9][4] = {
        {(1. - pktCal[0][0]) * (1. - pktCal[0][2]) / 4., (1. + pktCal[0][0]) * (1. - pktCal[0][2]) / 4., (1. + pktCal[0][0]) * (1. + pktCal[0][2]) / 4.,(1. - pktCal[0][0]) * (1. + pktCal[0][2]) / 4.  },
        {(1. - pktCal[1][0]) * (1. - pktCal[1][2]) / 4., (1. + pktCal[1][0]) * (1. - pktCal[1][2]) / 4., (1. + pktCal[1][0]) * (1. + pktCal[1][2]) / 4.,(1. - pktCal[1][0]) * (1. + pktCal[1][2]) / 4.  },
        {(1. - pktCal[2][0]) * (1. - pktCal[2][2]) / 4., (1. + pktCal[2][0]) * (1. - pktCal[2][2]) / 4., (1. + pktCal[2][0]) * (1. + pktCal[2][2]) / 4.,(1. - pktCal[2][0]) * (1. + pktCal[2][2]) / 4.  },
        {(1. - pktCal[3][0]) * (1. - pktCal[3][2]) / 4., (1. + pktCal[3][0]) * (1. - pktCal[3][2]) / 4., (1. + pktCal[3][0]) * (1. + pktCal[3][2]) / 4.,(1. - pktCal[3][0]) * (1. + pktCal[3][2]) / 4.  },
        {(1. - pktCal[4][0]) * (1. - pktCal[4][2]) / 4., (1. + pktCal[4][0]) * (1. - pktCal[4][2]) / 4., (1. + pktCal[4][0]) * (1. + pktCal[4][2]) / 4.,(1. - pktCal[4][0]) * (1. + pktCal[4][2]) / 4.  },
        {(1. - pktCal[5][0]) * (1. - pktCal[5][2]) / 4., (1. + pktCal[5][0]) * (1. - pktCal[5][2]) / 4., (1. + pktCal[5][0]) * (1. + pktCal[5][2]) / 4.,(1. - pktCal[5][0]) * (1. + pktCal[5][2]) / 4.  },
        {(1. - pktCal[6][0]) * (1. - pktCal[6][2]) / 4., (1. + pktCal[6][0]) * (1. - pktCal[6][2]) / 4., (1. + pktCal[6][0]) * (1. + pktCal[6][2]) / 4.,(1. - pktCal[6][0]) * (1. + pktCal[6][2]) / 4.  },
        {(1. - pktCal[7][0]) * (1. - pktCal[7][2]) / 4., (1. + pktCal[7][0]) * (1. - pktCal[7][2]) / 4., (1. + pktCal[7][0]) * (1. + pktCal[7][2]) / 4.,(1. - pktCal[7][0]) * (1. + pktCal[7][2]) / 4.  },
        {(1. - pktCal[8][0]) * (1. - pktCal[8][2]) / 4., (1. + pktCal[8][0]) * (1. - pktCal[8][2]) / 4., (1. + pktCal[8][0]) * (1. + pktCal[8][2]) / 4.,(1. - pktCal[8][0]) * (1. + pktCal[8][2]) / 4.  }
    };
};
//struct gauss { double macierz[nH * nH][nB * nB + 1] = {}; };
struct Sciany {
    double sciany[12][2] = {
        {-sqrt(3. / 5.), -1.}, //dol
        {0., -1.},
        {sqrt(3. / 5.), -1.},
        {1., -sqrt(3. / 5.)}, //prawo
        {1., 0.},
        {1., sqrt(3. / 5.)},
        {sqrt(3. / 5.), 1.}, //gora
        {0., 1.},
        {-sqrt(3. / 5.), 1.},
        {-1., sqrt(3. / 5.)}, //lewo
        {-1., 0.},
        {-1., -sqrt(3. / 5.)}
    };
    double funkcjeKsztaltu[12][4] = {
        {(1. - sciany[0][0]) * (1. - sciany[0][1]) / 4., (1. + sciany[0][0]) * (1. - sciany[0][1]) / 4., (1. + sciany[0][0]) * (1. + sciany[0][1]) / 4., (1. - sciany[0][0]) * (1. + sciany[0][1]) / 4.},
        {(1. - sciany[1][0]) * (1. - sciany[1][1]) / 4., (1. + sciany[1][0]) * (1. - sciany[1][1]) / 4., (1. + sciany[1][0]) * (1. + sciany[1][1]) / 4., (1. - sciany[1][0]) * (1. + sciany[1][1]) / 4.},
        {(1. - sciany[2][0]) * (1. - sciany[2][1]) / 4., (1. + sciany[2][0]) * (1. - sciany[2][1]) / 4., (1. + sciany[2][0]) * (1. + sciany[2][1]) / 4., (1. - sciany[2][0]) * (1. + sciany[2][1]) / 4.},

        {(1. - sciany[3][0]) * (1. - sciany[3][1]) / 4., (1. + sciany[3][0]) * (1. - sciany[3][1]) / 4., (1. + sciany[3][0]) * (1. + sciany[3][1]) / 4., (1. - sciany[3][0]) * (1. + sciany[3][1]) / 4.},
        {(1. - sciany[4][0]) * (1. - sciany[4][1]) / 4., (1. + sciany[4][0]) * (1. - sciany[4][1]) / 4., (1. + sciany[4][0]) * (1. + sciany[4][1]) / 4., (1. - sciany[4][0]) * (1. + sciany[4][1]) / 4.},
        {(1. - sciany[5][0]) * (1. - sciany[5][1]) / 4., (1. + sciany[5][0]) * (1. - sciany[5][1]) / 4., (1. + sciany[5][0]) * (1. + sciany[5][1]) / 4., (1. - sciany[5][0]) * (1. + sciany[5][1]) / 4.},

        {(1. - sciany[6][0]) * (1. - sciany[6][1]) / 4., (1. + sciany[6][0]) * (1. - sciany[6][1]) / 4., (1. + sciany[6][0]) * (1. + sciany[6][1]) / 4., (1. - sciany[6][0]) * (1. + sciany[6][1]) / 4.},
        {(1. - sciany[7][0]) * (1. - sciany[7][1]) / 4., (1. + sciany[7][0]) * (1. - sciany[7][1]) / 4., (1. + sciany[7][0]) * (1. + sciany[7][1]) / 4., (1. - sciany[7][0]) * (1. + sciany[7][1]) / 4.},
        {(1. - sciany[8][0]) * (1. - sciany[8][1]) / 4., (1. + sciany[8][0]) * (1. - sciany[8][1]) / 4., (1. + sciany[8][0]) * (1. + sciany[8][1]) / 4., (1. - sciany[8][0]) * (1. + sciany[8][1]) / 4.},
        
        {(1. - sciany[9][0]) * (1. - sciany[9][1]) / 4., (1. + sciany[9][0]) * (1. - sciany[9][1]) / 4., (1. + sciany[9][0]) * (1. + sciany[9][1]) / 4., (1. - sciany[9][0]) * (1. + sciany[9][1]) / 4.},
        {(1. - sciany[10][0]) * (1. - sciany[10][1]) / 4., (1. + sciany[10][0]) * (1. - sciany[10][1]) / 4., (1. + sciany[10][0]) * (1. + sciany[10][1]) / 4., (1. - sciany[10][0]) * (1. + sciany[10][1]) / 4.},
        {(1. - sciany[11][0]) * (1. - sciany[11][1]) / 4., (1. + sciany[11][0]) * (1. - sciany[11][1]) / 4., (1. + sciany[11][0]) * (1. + sciany[11][1]) / 4., (1. - sciany[11][0]) * (1. + sciany[11][1]) / 4.},

    };
};

struct macierz4x4 { double macierz[4][4] = {}; };

struct wektor { double wek[4] = {}; };

struct node {
    double x = 0., y = 0.;
    short int WB = 0;
};


struct Element {
    double data[3] = { -sqrt(3./5.), 0, sqrt(3./5.) };
    double NpoKSI[9][4] = {
        {(1. - data[0]) / (-4.), (1. - data[0]) / 4., (1. + data[0]) / 4., (1. + data[0]) / (-4.)},
        {(1. - data[0]) / (-4.), (1. - data[0]) / 4., (1. + data[0]) / 4., (1. + data[0]) / (-4.)},
        {(1. - data[0]) / (-4.), (1. - data[0]) / 4., (1. + data[0]) / 4., (1. + data[0]) / (-4.)},
        {(1. - data[1]) / (-4.), (1. - data[1]) / 4., (1. + data[1]) / 4., (1. + data[1]) / (-4.)},
        {(1. - data[1]) / (-4.), (1. - data[1]) / 4., (1. + data[1]) / 4., (1. + data[1]) / (-4.)},
        {(1. - data[1]) / (-4.), (1. - data[1]) / 4., (1. + data[1]) / 4., (1. + data[1]) / (-4.)},
        {(1. - data[2]) / (-4.), (1. - data[2]) / 4., (1. + data[2]) / 4., (1. + data[2]) / (-4.)},
        {(1. - data[2]) / (-4.), (1. - data[2]) / 4., (1. + data[2]) / 4., (1. + data[2]) / (-4.)},
        {(1. - data[2]) / (-4.), (1. - data[2]) / 4., (1. + data[2]) / 4., (1. + data[2]) / (-4.)}
    };
    double NpoETA[9][4] = {
        {(1. - data[0]) / (-4.), (1. + data[0]) / (-4.), (1. + data[0]) / 4., (1. - data[0]) / 4.},
        {(1. - data[1]) / (-4.), (1. + data[1]) / (-4.), (1. + data[1]) / 4., (1. - data[1]) / 4.},
        {(1. - data[2]) / (-4.), (1. + data[2]) / (-4.), (1. + data[2]) / 4., (1. - data[2]) / 4.},

        {(1. - data[0]) / (-4.), (1. + data[0]) / (-4.), (1. + data[0]) / 4., (1. - data[0]) / 4.},
        {(1. - data[1]) / (-4.), (1. + data[1]) / (-4.), (1. + data[1]) / 4., (1. - data[1]) / 4.},
        {(1. - data[2]) / (-4.), (1. + data[2]) / (-4.), (1. + data[2]) / 4., (1. - data[2]) / 4.},

        {(1. - data[0]) / (-4.), (1. + data[0]) / (-4.), (1. + data[0]) / 4., (1. - data[0]) / 4.},
        {(1. - data[1]) / (-4.), (1. + data[1]) / (-4.), (1. + data[1]) / 4., (1. - data[1]) / 4.},
        {(1. - data[2]) / (-4.), (1. + data[2]) / (-4.), (1. + data[2]) / 4., (1. - data[2]) / 4.}

    };
    /*double NpoKSI[9][4] = {
        {-0.443649167, 0.443649167, 0.056350833, -0.056350833},
        {-0.443649167, 0.443649167, 0.056350833, -0.056350833},
        {-0.443649167, 0.443649167, 0.056350833, -0.056350833},
        {-0.25,	0.25, 0.25, -0.25},
        {-0.25,	0.25, 0.25, -0.25},
        {-0.25,	0.25, 0.25, -0.25},
        {-0.056350833, 0.056350833,	0.443649167, -0.443649167},
        {-0.056350833, 0.056350833,	0.443649167, -0.443649167},
        {-0.056350833, 0.056350833,	0.443649167, -0.443649167}
    };
    double NpoETA[9][4] = {
        {-0.443649167, -0.056350833, 0.056350833, 0.443649167},
        {-0.25, -0.25, 0.25, 0.25},
        {-0.056350833, -0.443649167, 0.443649167, 0.056350833},

        {-0.443649167, -0.056350833, 0.056350833, 0.443649167},
        {-0.25, -0.25, 0.25, 0.25},
        {-0.056350833, -0.443649167, 0.443649167, 0.056350833},

        {-0.443649167, -0.056350833, 0.056350833, 0.443649167},
        {-0.25, -0.25, 0.25, 0.25},
        {-0.056350833, -0.443649167, 0.443649167, 0.056350833}

    };*/
};

struct Jakob {
    double J[2][2] = {};
    double Jinv[2][2] = {};
};

node* fillNodes(node t[]) {
    double dy = H / (nH - 1.);
    double dx = B / (nB - 1.);
    for (int i = 0; i < nB; i++) { // przechodze kolumnami
        for (int j = 0; j < nH; j++) {
            if (i == 0 || i == nB - 1 || j == 0 || j == nH - 1) t[j + i * nH].WB = 1;
            t[j + i * nH].x = (double)i * dx;
            t[j + i * nH].y = (double)j * dy;
            //cout << "<" << t[j + i * nH].x << ";" << t[j + i * nH].y << ";" << t[j + i * nH].WB << "> ";
        }
        //cout << endl;
    }
    return t;
}

macierz4x4 macierzHbc(Pkt p, macierz4x4 mac, Sciany s) {
    int licznik = 0;
    double L;
    double macierzPom[4][4] = {};
    for (int pkt = 0; pkt < 12; pkt += 3) {
        bool wejscDoPetli = false;
        switch (licznik) {
            //SPRAWDZAM CZY DANA SCIANA MOZE WEJSC DO PETLI
        case 0: if (p.pkt[0][2] == 1 && p.pkt[1][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[1][0] - p.pkt[0][0], 2) + pow(p.pkt[1][1] - p.pkt[0][1], 2)); break;
        case 1: if (p.pkt[1][2] == 1 && p.pkt[2][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[2][0] - p.pkt[1][0], 2) + pow(p.pkt[2][1] - p.pkt[1][1], 2)); break;
        case 2: if (p.pkt[2][2] == 1 && p.pkt[3][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[3][0] - p.pkt[2][0], 2) + pow(p.pkt[3][1] - p.pkt[2][1], 2)); break;
        case 3: if (p.pkt[3][2] == 1 && p.pkt[0][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[0][0] - p.pkt[3][0], 2) + pow(p.pkt[0][1] - p.pkt[3][1], 2)); break;
        }
        if (wejscDoPetli) {
            double detJ = L / 2.;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    //LICZE MACIERZ * MACIERZ TRANSP
                    macierzPom[i][j]  = 5. / 9. * s.funkcjeKsztaltu[pkt][j] * s.funkcjeKsztaltu[pkt][i] * alpha;
                    macierzPom[i][j]  += 8. / 9. * s.funkcjeKsztaltu[pkt +1][j] * s.funkcjeKsztaltu[pkt +1][i] * alpha;
                    macierzPom[i][j]  += 5. / 9. * s.funkcjeKsztaltu[pkt +2][j] * s.funkcjeKsztaltu[pkt +2][i] * alpha;

                    //OBLICZAM MACIERZE HBCi I DODAJE DO MACIERZY KONCOWEJ
                    mac.macierz[i][j] += macierzPom[i][j] * detJ;
                }
            }
        }
        licznik++;
    }
    return mac;
}

macierz4x4 macierzH(macierz4x4 mac, Pkt p, Element e, Jakob j) {
    double dNpodX[9][4] = {};
    double dNpodY[9][4] = {};
    double NpoX[9][4] = {};
    double NpoY[9][4] = {};
    double det;
    for (int pkt = 0; pkt < 9; pkt++) {
        for (int i = 0; i < 4; i++) {
            //dxpodksi
            j.J[0][0] += e.NpoKSI[pkt][i] * p.pkt[i][0];
            //dxpodeta
            j.J[1][0] += e.NpoETA[pkt][i] * p.pkt[i][0];
            //dypodksi
            j.J[0][1] += e.NpoKSI[pkt][i] * p.pkt[i][1];
            //dypodeta
            j.J[1][1] += e.NpoETA[pkt][i] * p.pkt[i][1];
        }
        if (j.J[0][0] < DOKLADNOSC) j.J[0][0] = 0.;
        if (j.J[1][0] < DOKLADNOSC) j.J[1][0] = 0.;
        if (j.J[0][1] < DOKLADNOSC) j.J[0][1] = 0.;
        if (j.J[1][1] < DOKLADNOSC) j.J[1][1] = 0.;

        det = j.J[0][0] * j.J[1][1] - j.J[0][1] * j.J[1][0];
        j.Jinv[0][0] = j.J[1][1] / det; // * 1/det
        j.Jinv[0][1] = j.J[0][1] / -det;
        j.Jinv[1][0] = j.J[1][0] / -det; 
        j.Jinv[1][1] = j.J[0][0] / det;

        //TWORZE MACIERZE dNi po dX i dY
        for (int x = 0; x < 9; x++) {
            for (int y = 0; y < 4; y++) {
                dNpodX[x][y] = j.Jinv[0][0] * e.NpoKSI[x][y] + j.Jinv[1][0] * e.NpoETA[x][y];
                dNpodY[x][y] = j.Jinv[0][1] * e.NpoKSI[x][y] + j.Jinv[1][1] * e.NpoETA[x][y];
            }
        }
        //mnoze macierz * macierz^T
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                NpoX[x][y] = dNpodX[pkt][x] * dNpodX[pkt][y];
                NpoY[x][y] = dNpodY[pkt][x] * dNpodY[pkt][y];
                //licze macierz H
                mac.macierz[x][y] += kodt * p.pktCal[pkt][1] * p.pktCal[pkt][3] * (NpoX[x][y] + NpoY[x][y]) * det;               
            }
        }
        j.J[0][0] = 0.;
        j.J[0][1] = 0.;
        j.J[1][0] = 0.;
        j.J[1][1] = 0.;
    }
    return mac;
}

wektor wektP(Pkt p, wektor w, Sciany s) {
    int licznik = 0;
    double L;
    double wektorPom[4] = {};
    double detJ;
    for (int pkt = 0; pkt < 12; pkt += 3) {
        bool wejscDoPetli = false;
        switch (licznik) {
            //SPRAWDZAM CZY DANA SCIANA MOZE WEJSC DO PETLI
        case 0: if (p.pkt[0][2] == 1 && p.pkt[1][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[1][0] - p.pkt[0][0], 2) + pow(p.pkt[1][1] - p.pkt[0][1], 2)); break;
        case 1: if (p.pkt[1][2] == 1 && p.pkt[2][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[2][0] - p.pkt[1][0], 2) + pow(p.pkt[2][1] - p.pkt[1][1], 2)); break;
        case 2: if (p.pkt[2][2] == 1 && p.pkt[3][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[3][0] - p.pkt[2][0], 2) + pow(p.pkt[3][1] - p.pkt[2][1], 2)); break;
        case 3: if (p.pkt[3][2] == 1 && p.pkt[0][2] == 1) wejscDoPetli = true; else wejscDoPetli = false;
            L = sqrt(pow(p.pkt[0][0] - p.pkt[3][0], 2) + pow(p.pkt[0][1] - p.pkt[3][1], 2)); break;
        }
        if (wejscDoPetli) {
            detJ = L / 2.;
            for (int i = 0; i < 4; i++) {
                //LICZE waga * wektor * temperatura otoczenia - wzor z pdf
                wektorPom[i] = 5. / 9. * s.funkcjeKsztaltu[pkt][i] * tempOt;
                wektorPom[i] += 8. / 9. * s.funkcjeKsztaltu[pkt + 1][i] * tempOt;
                wektorPom[i] += 5. / 9. * s.funkcjeKsztaltu[pkt + 2][i] * tempOt;
                //licze macierz P lokalna
                //sumuje obie macierze, mnoze przez wsp (WZOR) - dodaje wszystko do ostatecznej macierzy
                w.wek[i] += alpha * wektorPom[i] * detJ;
            }
        }
        licznik++;
    }
    return w;
}

macierz4x4 macierzC(macierz4x4 mac, Pkt p, Element e, Jakob j) {
    double det;
    for (int pkt = 0; pkt < 9; pkt++) {
        for (int i = 0; i < 4; i++) {
            //dxpodksi
            j.J[0][0] += e.NpoKSI[pkt][i] * p.pkt[i][0];
            //dxpodeta
            j.J[1][0] += e.NpoETA[pkt][i] * p.pkt[i][0];
            //dypodksi
            j.J[0][1] += e.NpoKSI[pkt][i] * p.pkt[i][1];
            //dypodeta
            j.J[1][1] += e.NpoETA[pkt][i] * p.pkt[i][1];
        }
        if (j.J[0][0] < DOKLADNOSC) j.J[0][0] = 0.;
        if (j.J[1][0] < DOKLADNOSC) j.J[1][0] = 0.;
        if (j.J[0][1] < DOKLADNOSC) j.J[0][1] = 0.;
        if (j.J[1][1] < DOKLADNOSC) j.J[1][1] = 0.;

        det = j.J[0][0] * j.J[1][1] - j.J[0][1] * j.J[1][0];
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                //licze macierz C
                mac.macierz[x][y] += c * ro * p.pktCal[pkt][1] * p.pktCal[pkt][3] * p.funkcjeKsztaltu[pkt][x] * p.funkcjeKsztaltu[pkt][y] * det;
            }
        }
        j.J[0][0] = 0.;
        j.J[0][1] = 0.;
        j.J[1][0] = 0.;
        j.J[1][1] = 0.;
    }
    return mac;
}

void Gauss(double* g[], double* temp) {
    int i, j, k;
    double m, s;
    /*double AB[nH * nH][nB * nB + 1] = {};
    double T1[nB * nB] = {};

    for (int i = 0; i < (nH * nH); i++) {
        for (int j = 0; j < (nB * nB + 1); j++) {
            AB[i][j] = g.macierz[i][j];
        }
    }*/
    double** AB = new double* [nH * nH];// [nH * nH] [nB * nB + 1] = {};
    for (int i = 0; i < nH * nH; i++) {
        AB[i] = new double[nH * nH + 1];
    }
    //double T1[nB * nB] = {};
    double* T1 = new double[nH * nH];

    for (int i = 0; i < (nH * nH); i++) {
        for (int j = 0; j < (nB * nB + 1); j++) {
            //AB[i][j] = g.macierz[i][j];
            AB[i][j] = g[i][j];
        }
        T1[i] = 0.;
    }
    // eliminacja wspó³czynników
    for (i = 0; i < (nH * nH) - 1; i++)
    {
        for (j = i + 1; j < (nB * nB); j++)
        {
            if (fabs(AB[i][i]) < DOKLADNOSC) {
                std::cout << "\n\nDzielenie przez 0.\n\n";
                return;
            }
            m = -AB[j][i] / AB[i][i];
            for (k = i + 1; k <= (nH * nH); k++)
                AB[j][k] += m * AB[i][k];
        }
    }

    // wyliczanie niewiadomych
    for (i = (nH * nH) - 1; i >= 0; i--)
    {
        s = AB[i][(nH * nH)];
        for (j = (nH * nH) - 1; j >= i + 1; j--)
            s -= AB[i][j] * T1[j];
        if (fabs(AB[i][i]) < DOKLADNOSC) {
            std::cout << "\n\nDzielenie przez 0.\n\n";
            return;
        }
        T1[i] = s / AB[i][i];
    }

    double minT = 1.e20;
    double maxT = 0.;
    for (i = 0; i < nH * nB; i++) {
        temp[i] = T1[i];
        //cout << setw(10) << T1[i] << endl;
        if (T1[i] > maxT) maxT = T1[i];
        if (T1[i] < minT) minT = T1[i];
    }
    std::cout << "MinT: " << minT;
    std::cout << "\nMaxT: " << maxT << endl;
}

int main() {
    node pom[nH * nB] = {};
    node* nodes = fillNodes(pom);
    Pkt pkt;
    macierz4x4 MAC;
    wektor WEK;
    Element el;
    Jakob j;
    Sciany s;
    //gauss g;
    int ID[4] = {};
    //double wektorTemperatur[nH * nB] = {};
    double* wektorTemperatur = new double[nH * nB];

    // wektor temperatury poczatkowej t0
    for (int i = 0; i < nB * nB; i++) {
        wektorTemperatur[i] = t0;
    }
    /*double agregacjaH[nH * nB][nH * nB] = {};
    double agregacjaC[nH * nB][nH * nB] = {};
    double macierzCrazywektorT0[nH * nB] = {};
    double wekP[nH * nB] = {};*/

    double** agregacjaH = new double* [nH * nB];
    double** agregacjaC = new double* [nH * nB];
    for (int i = 0; i < nH * nB; i++) {
        agregacjaH[i] = new double[nH * nB];
        agregacjaC[i] = new double[nH * nB];
    }

    double* wekP = new double[nH * nB];
    double* macierzCrazywektorT0 = new double[nH * nB];
    for (int i = 0; i < nH * nB; i++) {
        for (int j = 0; j < nH * nB; j++) {
            agregacjaC[i][j] = 0.;
            agregacjaH[i][j] = 0.;
        }
        wekP[i] = 0.;
        macierzCrazywektorT0[i] = 0.;
    }

    double** macierzGauss = new double* [nH * nH];
    for (int i = 0; i < nH * nH; i++)
        macierzGauss[i] = new double[nH * nH + 1];
    for (int i = 0; i < nH * nH; i++) {
        for (int j = 0; j < nB * nB + 1; j++) {
            macierzGauss[i][j] = 0.;
        }
    }

    wektor wektorP;
    macierz4x4 macC;
    macierz4x4 macH;
    macierz4x4 macHBC;

    for (int czas = 0; czas < t; czas += dt) {
        cout << "Czas " << czas+dt << endl;
        int temp = 0, temp2 = 1; //zmienne sterujace petla po elementach
        for (int i = 0; i < (nH - 1) * (nB - 1); i++) {
            //Licze wartosci ID wezlow ktore tworza element
            ID[0] = i + temp;
            ID[1] = i + temp + nH;
            ID[2] = i + temp + nH + 1;
            ID[3] = i + temp + 1;
            //Wybieram 4 punkty tworzace element
            pkt.pkt[0][0] = nodes[ID[0]].x; pkt.pkt[0][1] = nodes[ID[0]].y; pkt.pkt[0][2] = nodes[ID[0]].WB;//x1y1
            pkt.pkt[1][0] = nodes[ID[1]].x; pkt.pkt[1][1] = nodes[ID[1]].y; pkt.pkt[1][2] = nodes[ID[1]].WB;//x2y2
            pkt.pkt[2][0] = nodes[ID[2]].x; pkt.pkt[2][1] = nodes[ID[2]].y; pkt.pkt[2][2] = nodes[ID[2]].WB;//x3y3
            pkt.pkt[3][0] = nodes[ID[3]].x; pkt.pkt[3][1] = nodes[ID[3]].y; pkt.pkt[3][2] = nodes[ID[3]].WB;//x4y4

            macH = macierzH(MAC, pkt, el, j);
            macHBC = macierzHbc(pkt, MAC, s);
            macC = macierzC(MAC, pkt, el, j);
            wektorP = wektP(pkt, WEK, s);

            for (int x = 0; x < 4; x++) {
                for (int y = 0; y < 4; y++) {
                    // [H]+[Hbc]+[C]
                    agregacjaH[ID[x]][ID[y]] += macH.macierz[x][y] + macHBC.macierz[x][y] + macC.macierz[x][y] / dt;
                    agregacjaC[ID[x]][ID[y]] += macC.macierz[x][y] / dt;
                }
                // {P}
                wekP[ID[x]] += wektorP.wek[x];
            }
            temp2++;
            if (temp2 == nH) {
                temp2 = 1;
                temp++;
            }
        }

        for (int i = 0; i < nH * nB; i++) {
            for (int j = 0; j < nH * nB; j++) {
                // {[C]/dt} * {T0}
                macierzCrazywektorT0[i] += agregacjaC[i][j] * wektorTemperatur[j];
            }
            // {P} + {[C]/dt} * {T0}
            wekP[i] += macierzCrazywektorT0[i];
        }


        for (int i = 0; i < nH * nH; i++) {
            for (int j = 0; j < nB * nB + 1; j++) {
                //if (j != nB*nB) g.macierz[i][j] = agregacjaH[i][j];
                if (j != nB * nB) macierzGauss[i][j] = agregacjaH[i][j];
                //else g.macierz[i][j] = wekP[i];
                else macierzGauss[i][j] = wekP[i];
            }
        }
        /*for (int i = 0; i < nB * nH; i++) {
            for (int j = 0; j < nH * nB; j++) {
                cout << setw(10);
                cout << agregacjaH[i][j];
            }
            cout << " |  " << wekP[i] << endl;
            cout << endl;
        }*/
        Gauss(macierzGauss, wektorTemperatur);
        //wyzerowanie macierzy
        for (int i = 0; i < nH * nH; i++) {
            for (int j = 0; j < nB * nB; j++) {
                agregacjaH[i][j] = 0.;
                agregacjaC[i][j] = 0.;
            }
            macierzCrazywektorT0[i] = 0.;
            wekP[i] = 0.;
        }
    }
    delete[] wektorTemperatur;
    delete[] wekP;
    delete[] macierzGauss;
    delete[] agregacjaC;
    delete[] agregacjaH;
    delete[] macierzCrazywektorT0;

}
