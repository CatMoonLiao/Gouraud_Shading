#include <iostream>
#include <GL/glut.h>
#include <string>
#include <fstream>
#include <math.h>
#define PI 3.141592
int WIN_WID = 500;
int WIN_HEI = 500;
using namespace std;
FILE *file;
struct point {
    float x;
    float y;
    float z;
    float w;
};

struct color {
    float R;
    float G;
    float B;
};

struct plane {
    int vertice_c;
    int vertices[4];
    point n;
};

typedef struct object *obj_ptr;
struct object {
    bool exist;
    string fname;
    int pointc;
    int planec;
    float obj_point[4][2000] = { 0 };
    plane obj_plane[3800] = { 0 };
    point point_normal[2000] = { 0 };
    color point_color[2000] = { 0 };
};


typedef struct Data *data_ptr;
struct Data {
    object obj;
    data_ptr link;
    color O;
    float Ks;
    float Kd;
    int N;
};
data_ptr datalist = 0, before = 0;

struct light {
    color lightcolor;
    point location;
};
light lightlist[10] = { 0 };
int lightc = 0;

float **TM;
float **TEMP;
float **EM;
float **PM;
point wmax, wmin;
point vmin, vmax, WIN_PC, pc;
float AR = 1;

point node[2] = { 0 };
color nodec[2] = { 0 };
point E, COI;
object skull = { 0 }, teapot = { 0 }, cube = { 0 }, grid = { 0 };
float tilt;
float PMh, PMy, PMtheta;
bool first_obj = true;
color bgcolor, kaia;
float ZB[600][600];
color CB[600][600];

void dot(int x, int y) {
    glBegin(GL_POINTS);
    glVertex2i(x, WIN_HEI - y);
    glEnd();
}

//sqare matrix mul
void mul(float **A, float **B) {

    float ANS[4][4] = { 0 };
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                ANS[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    //copy
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            B[i][j] = ANS[i][j];
}

//initial matrix to unit matrix
void reset(float **A) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            A[i][j] = 0;
            if (i == j)
                A[i][j] = 1;
        }
    }
}

void Initial(void)//初始化函数 
{
    wmax.x = 1;
    wmax.y = 1;
    wmin.x = -1;
    wmin.y = -1;
    //initial TM
    TM = new float *[4];
    for (int i = 0; i < 4; i++) {
        TM[i] = new float[4];
    }
    reset(TM);
    //initial TEMP
    TEMP = new float *[4];
    for (int i = 0; i < 4; i++) {
        TEMP[i] = new float[4];
    }
    reset(TEMP);
    //initial EM
    EM = new float *[4];
    for (int i = 0; i < 4; i++) {
        EM[i] = new float[4];
    }
    reset(EM);
    //initial PM
    PM = new float *[4];
    for (int i = 0; i < 4; i++) {
        PM[i] = new float[4];
    }
    reset(PM);

    glColor3f(0, 0, 0);
    glPointSize(1.0);
    glMatrixMode(GL_PROJECTION);
    glutInitWindowSize(WIN_WID, WIN_HEI);//设定窗口的大小
    gluOrtho2D(0, WIN_WID, WIN_HEI, 0);
    glClearColor(0, 0, 0, 0);//背景顏色
    glClear(GL_COLOR_BUFFER_BIT);//用背景色清空畫布
    glFinish();
}

void myDisplay(void)//显示回调函数
{
    glFlush();//清空OpenGL命令缓冲区，强制执行命令缓冲区中所有OpenGL函数
}

//將座標系統轉為pixel 並4捨5入
point coord2pixel(point a, float x, float y, point c) {
    a.x *= x;
    a.y *= y;
    a.x += c.x;
    a.y += c.y;
    a.x = floor(a.x + 0.5);
    a.y = floor(a.y + 0.5);

    return a;
}

point Normalize(point a) {
    float length = a.x*a.x + a.y*a.y + a.z*a.z;
    length = sqrt(length);
    a.x /= length;
    a.y /= length;
    a.z /= length;
    return a;
}

void PerspectiveDivide(point *a) {
    a->x /= a->w;
    a->y /= a->w;
    a->z /= a->w;
    a->w /= a->w;
}

point Cross(point a, point b) {
    point ans;
    ans.x = a.y*b.z - a.z*b.y;
    ans.y = a.z*b.x - a.x*b.z;
    ans.z = a.x*b.y - a.y*b.x;
    return ans;
}

float Dot(point a, point b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

float max(float a, float b) {
    if (a > b)
        return a;
    else
        return b;
}

float min(float a, float b) {
    if (a < b)
        return a;
    else
        return b;
}

void mul_obj(float **A, obj_ptr o) {
    float ans[4][2000] = { 0 };
    int col = o->pointc;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < col; j++) {
            for (int k = 0; k < 4; k++) {
                ans[i][j] += A[i][k] * o->obj_point[k][j];
            }
        }
    }
    //copy
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < col; j++)
            o->obj_point[i][j] = ans[i][j];
}

//position是用來算的點 normal是法向量
color illuminate(color O, float Kd, float Ks, int N, point position, point normal) {
    color ans = kaia;
    normal = Normalize(normal);
    //ambient
    ans.R *= O.R;
    ans.G *= O.G;
    ans.B *= O.B;
    //diffuse
    color diff = { 0 }, spec = { 0 };
    for (int i = 0; i < lightc; i++) {
        light templight = lightlist[i];
        point L;
        L.x = templight.location.x - position.x;
        L.y = templight.location.y - position.y;
        L.z = templight.location.z - position.z;
        L = Normalize(L);
        //diffuse
        float dotv = Dot(L, normal);
        if (dotv < 0)
            dotv = 0;

        diff.R += Kd * templight.lightcolor.R*dotv;
        diff.G += Kd * templight.lightcolor.G*dotv;
        diff.B += Kd * templight.lightcolor.B*dotv;
        //specular
        //角平分線
        point H, V;
        V.x = E.x - position.x;
        V.y = E.y - position.y;
        V.z = E.z - position.z;
        V = Normalize(V);

        H.x = V.x + L.x;
        H.y = V.y + L.y;
        H.z = V.z + L.z;
        H = Normalize(H);
        dotv = Dot(H, normal);
        dotv = pow(dotv, N);
        if (dotv < 0)
            dotv = 0;
        spec.R += Ks * templight.lightcolor.R*dotv;
        spec.G += Ks * templight.lightcolor.G*dotv;
        spec.B += Ks * templight.lightcolor.B*dotv;
    }
    diff.R *= O.R;
    diff.G *= O.G;
    diff.B *= O.B;
    ans.R += diff.R + spec.R;
    ans.G += diff.G + spec.G;
    ans.B += diff.B + spec.B;
    //clamp

    return ans;
}

bool inPolygon(int planei, point target, point  p[], int vertice_c) {
    bool clock = true;
    point v1, v2;
    v1.x = p[1].x - p[0].x;
    v1.y = p[1].y - p[0].y;
    v1.z = p[1].z - p[0].z;
    v2.x = p[2].x - p[1].x;
    v2.y = p[2].y - p[1].y;
    v2.z = p[2].z - p[1].z;

    float f = Cross(v1, v2).z;
    if (f > 0)
        clock = false;
    for (int i = 0; i < vertice_c; i++) {
        int next = i + 1;
        if (i + 1 == vertice_c)
            next = 0;

        int t = (target.x - p[i].x)*(p[next].y - p[i].y) - (p[next].x - p[i].x)*(target.y - p[i].y);
        if (clock) {
            if (t < 0)
                return false;
        }
        else {
            if (t > 0)
                return false;
        }
    }
    return true;
}

float findZ(point target, point normal, point p) {
    //a(x-x0)+b(y-y0)+c(z-z0)=0
    return p.z - ((target.x - p.x)*normal.x + (target.y - p.y)*normal.y) / normal.z;
}

int onLine(int vc, int targetY, color pc[], point p[]) {
    int c = 0;
    for (int i = 0; i < vc; i++) {

        int next = i + 1;
        if (next == vc)
            next = 0;
        //cout << pc[i].R << " " << pc[i].G << " " << pc[i].B << "||" << pc[next].R << " " << pc[next].G << " " << pc[next].B << endl;
        //system("pause");
        if (p[next].y - p[i].y == 0) {
            if (p[next].y != targetY)
                continue;
            else {
                nodec[c] = pc[i];
                node[c++] = p[i];
                nodec[c] = pc[next];
                node[c++] = p[next];
                return c;
            }
        }
        float x = (p[next].x - p[i].x)*(targetY - p[i].y) / (p[next].y - p[i].y) + p[i].x;
        float t = min(p[next].y, p[i].y);
        if (t > targetY)
            continue;
        t = max(p[next].y, p[i].y);
        if (t < targetY)
            continue;

        node[c].x = floor(x + 0.5);
        node[c].y = targetY;
        nodec[c].R = pc[i].R + (pc[next].R - pc[i].R)*(targetY - p[i].y) / (p[next].y - p[i].y);
        nodec[c].G = pc[i].G + (pc[next].G - pc[i].G)*(targetY - p[i].y) / (p[next].y - p[i].y);
        nodec[c].B = pc[i].B + (pc[next].B - pc[i].B)*(targetY - p[i].y) / (p[next].y - p[i].y);
        c++;
        if (c > 1 && node[0].x == node[1].x)
            c--;
        if (c == 2)
            return 2;
    }
    return c;
}

void draw(object o) {
    for (int k = 0; k < o.planec; k++) {
        /*取點的值*/
        point pvalue[4] = { 0 };
        color pcolor[4] = { 0 };
        for (int j = 0; j < o.obj_plane[k].vertice_c; j++) {
            pvalue[j].x = o.obj_point[0][o.obj_plane[k].vertices[j] - 1];
            pvalue[j].y = o.obj_point[1][o.obj_plane[k].vertices[j] - 1];
            pvalue[j].z = o.obj_point[2][o.obj_plane[k].vertices[j] - 1];
            pvalue[j].w = o.obj_point[3][o.obj_plane[k].vertices[j] - 1];
            PerspectiveDivide(&pvalue[j]);
            pvalue[j] = coord2pixel(pvalue[j], (vmax.x - vmin.x) / 2, (vmax.y - vmin.y) / 2, pc);
            pcolor[j] = o.point_color[o.obj_plane[k].vertices[j] - 1];
        }
        point v1, v2;
        v1.x = pvalue[1].x - pvalue[0].x;
        v1.y = pvalue[1].y - pvalue[0].y;
        v1.z = pvalue[1].z - pvalue[0].z;
        v2.x = pvalue[2].x - pvalue[1].x;
        v2.y = pvalue[2].y - pvalue[1].y;
        v2.z = pvalue[2].z - pvalue[1].z;
        point normal = Cross(v2, v1);
        /*求平面的xy最大值*/
        point maxp, minp;
        maxp.x = 0;
        maxp.y = 0;
        minp.x = 1000;
        minp.y = 1000;
        for (int i = 0; i < o.obj_plane[k].vertice_c; i++) {
            maxp.x = max(pvalue[i].x, maxp.x);
            maxp.y = max(pvalue[i].y, maxp.y);
            minp.x = min(pvalue[i].x, minp.x);
            minp.y = min(pvalue[i].y, minp.y);
        }

        /*算橫排*/
        for (int i = maxp.y; i >= minp.y; i--) {
            node[0].x = 0;
            node[0].y = 0;
            node[1].x = 0;
            node[1].y = 0;
            nodec[0].R = 0;
            nodec[0].G = 0;
            nodec[0].B = 0;
            nodec[1].R = 0;
            nodec[1].G = 0;
            nodec[1].B = 0;

            int ncount = onLine(o.obj_plane[k].vertice_c, i, pcolor, pvalue);
            //按比例算整排
            if (ncount == 1) {
                float z = findZ(node[0], normal, pvalue[0]);
                int x = node[0].x, y = node[0].y;
                if (z < ZB[y][x]) {
                    ZB[y][x] = z;
                    CB[y][x] = nodec[0];
                }
            }
            else {
                if (node[0].x > node[1].x) {
                    swap(node[0], node[1]);
                    swap(nodec[0], nodec[1]);
                }
                int dx = node[1].x - node[0].x;
                for (int j = node[0].x; j <= node[1].x; j++) {
                    point target;
                    target.x = j;
                    target.y = i;
                    float z = findZ(target, normal, pvalue[1]);
                    if (z < ZB[i][j]) {
                        ZB[i][j] = z;
                        color temp = { 0 };
                        temp.R = nodec[0].R + (nodec[1].R - nodec[0].R)*(j - node[0].x) / dx;
                        temp.G = nodec[0].G + (nodec[1].G - nodec[0].G)*(j - node[0].x) / dx;
                        temp.B = nodec[0].B + (nodec[1].B - nodec[0].B)*(j - node[0].x) / dx;
                        CB[i][j] = temp;
                    }
                }
            }
        }

    }

}

void printM(float **A) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

void culallnormal(data_ptr dp) {
    //every plane normal
    for (int i = 0; i < dp->obj.planec; i++) {
        point v1, v2;
        int i1 = dp->obj.obj_plane[i].vertices[0] - 1;
        int i2 = dp->obj.obj_plane[i].vertices[1] - 1;
        int i3 = dp->obj.obj_plane[i].vertices[2] - 1;
        v1.x = dp->obj.obj_point[0][i2] - dp->obj.obj_point[0][i1];
        v1.y = dp->obj.obj_point[1][i2] - dp->obj.obj_point[1][i1];
        v1.z = dp->obj.obj_point[2][i2] - dp->obj.obj_point[2][i1];
        v2.x = dp->obj.obj_point[0][i3] - dp->obj.obj_point[0][i2];
        v2.y = dp->obj.obj_point[1][i3] - dp->obj.obj_point[1][i2];
        v2.z = dp->obj.obj_point[2][i3] - dp->obj.obj_point[2][i2];
        dp->obj.obj_plane[i].n = Cross(v2, v1);
    }

    //every vertices normal
    int nc[2000] = { 0 };
    point vn[2000] = { 0 };
    for (int i = 0; i < dp->obj.planec; i++) {
        for (int j = 0; j < dp->obj.obj_plane[i].vertice_c; j++) {
            int a = dp->obj.obj_plane[i].vertices[j] - 1;
            nc[a]++;
            vn[a].x += dp->obj.obj_plane[i].n.x;
            vn[a].y += dp->obj.obj_plane[i].n.y;
            vn[a].z += dp->obj.obj_plane[i].n.z;
        }
    }
    //save vertices normal value
    for (int i = 0; i < dp->obj.pointc; i++) {
        dp->obj.point_normal[i].x = vn[i].x / nc[i];
        dp->obj.point_normal[i].y = vn[i].y / nc[i];
        dp->obj.point_normal[i].z = vn[i].z / nc[i];
    }
}

//***calculate vertices color***//
void culvertexcolor(data_ptr dp) {
    for (int i = 0; i < dp->obj.pointc; i++) {
        point v = { 0 };
        v.x = dp->obj.obj_point[0][i];
        v.y = dp->obj.obj_point[1][i];
        v.z = dp->obj.obj_point[2][i];
        v.w = dp->obj.obj_point[3][i];
        dp->obj.point_color[i] = illuminate(dp->O, dp->Kd, dp->Ks, dp->N, v, dp->obj.point_normal[i]);
    }
}

void View() {
    reset(TEMP);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            TEMP[i][j] = EM[i][j];
    mul(PM, TEMP);

    //initial z-buffer
    for (int i = 0; i < 600; i++) {
        for (int j = 0; j < 600; j++) {
            //view
            CB[i][j] = { 0 };
            ZB[i][j] = 0;
            if (i <= vmax.y&&i >= vmin.y) {
                if (j <= vmax.x&&j >= vmin.x) {
                    CB[i][j] = bgcolor;
                    ZB[i][j] = 1;
                }
            }
        }
    }
    data_ptr go = datalist;
    while (go) {
        culallnormal(go);
        culvertexcolor(go);
        object o = go->obj;
        mul_obj(TEMP, &o);
        draw(o);
        go = go->link;
    }
    for (int i = 0; i < WIN_HEI; i++) {
        for (int j = 0; j < WIN_WID; j++) {
            glColor3f(CB[i][j].R, CB[i][j].G, CB[i][j].B);
            dot(j, i);
        }
    }
    glFinish();
}

obj_ptr readASC(string fn) {
    FILE *f;
    errno_t err = fopen_s(&f, fn.c_str(), "r");
    if (err != 0) {
        cout << "Can't open this asc file." << endl;
        return 0;
    }

    object temp = { 0 };
    temp.planec = 100;
    temp.fname = fn;
    char c = 0;
    string in = "";
    int linec = 0;
    while (c = fgetc(f)) {
        if (c != '\n' &&c != EOF)
            in += c;
        else {
            if (in.length() == 0)
                continue;//跳過空白行
            linec++;
            if (linec == 1) {
                int sum = 0;
                bool p = true;
                for (int i = 0; i < in.length(); i++) {
                    if (in[i] >= '0' && in[i] <= '9') {
                        sum *= 10;
                        sum += in[i] - '0';
                    }
                    else {
                        if (p) {
                            temp.pointc = sum;
                            sum = 0;
                            p = false;
                        }
                        else {
                            if (in[i - 1] != ' ')
                                break;
                        }
                    }
                }
                temp.planec = sum;
            }
            else {
                //input is point
                if (linec <= temp.pointc + 1) {
                    int index = linec - 2;
                    int count = 0;
                    float sum = 0, fordot = 1;
                    bool dot = false;
                    for (int i = 0; i < in.length(); i++) {

                        if (in[i] >= '0' && in[i] <= '9') {
                            if (dot)
                                fordot *= 10;
                            sum *= 10;
                            sum += in[i] - '0';
                            continue;
                        }
                        else if (in[i] == '-') {
                            fordot *= -1;
                            continue;
                        }
                        else if (in[i] == '.') {
                            dot = true;
                            continue;
                        }
                        else {//space
                            if (i > 0 && in[i - 1] != ' ') {//非連續空白
                                sum /= fordot;
                                switch (count) {
                                case 0:
                                    temp.obj_point[count++][index] = sum;
                                    break;
                                case 1:
                                    temp.obj_point[count++][index] = sum;
                                    break;
                                case 2:
                                    temp.obj_point[count++][index] = sum;
                                    break;
                                }
                                dot = false;
                                sum = 0;
                                fordot = 1;
                            }
                        }
                        if (count == 2 && i == in.length() - 1)
                            break;
                        else if (count == 3)
                            break;
                    }
                    temp.obj_point[count][index] = sum / fordot;
                    temp.obj_point[3][index] = 1;
                }
                //input is plane
                else {

                    int index = linec - temp.pointc - 2;
                    int count = 0, sum = 0;
                    for (int i = 0; i < in.length(); i++) {
                        if (i == 0) {
                            switch (in[i]) {
                            case '3':
                                temp.obj_plane[index].vertice_c = 3;
                                break;
                            case '4':
                                temp.obj_plane[index].vertice_c = 4;
                                break;
                            }
                            continue;
                        }

                        if (in[i] <= '9' &&in[i] >= '0') {
                            sum *= 10;
                            sum += in[i] - '0';
                            if (i == in.length() - 1) {
                                temp.obj_plane[index].vertices[count] = sum;
                                break;
                            }
                        }
                        else if (i > 1 && in[i] == ' ' && in[i - 1] != ' ') {
                            temp.obj_plane[index].vertices[count++] = sum;
                            sum = 0;
                        }
                    }
                }
            }
            in = "";
        }
        if (c == EOF || linec > temp.planec + temp.pointc) {
            cout << "read asc file end." << endl;
            break;
        }
    }
    fclose(f);
    temp.exist = true;
    if (fn == "teapot.asc") {
        teapot = temp;
        return &teapot;
    }
    else {
        grid = temp;
        return &grid;
    }
}

void Translate(point a, float **A) {
    reset(TEMP);
    TEMP[0][3] = a.x;
    TEMP[1][3] = a.y;
    TEMP[2][3] = a.z;
    mul(TEMP, A);
}

void Scale(point a) {
    reset(TEMP);
    TEMP[0][0] = a.x;
    TEMP[1][1] = a.y;
    TEMP[2][2] = a.z;
    mul(TEMP, TM);
}

void Rotate(point a) {
    reset(TEMP);
    a.x *= PI / 180;
    TEMP[1][1] = (float)cos(a.x);
    TEMP[2][1] = (float)sin(a.x);
    TEMP[1][2] = (float)sin(a.x)*(-1);
    TEMP[2][2] = (float)cos(a.x);
    mul(TEMP, TM);

    reset(TEMP);
    a.y *= PI / 180;
    TEMP[0][0] = (float)cos(a.y);
    TEMP[2][0] = (float)sin(a.y)*(-1);
    TEMP[0][2] = (float)sin(a.y);
    TEMP[2][2] = (float)cos(a.y);
    mul(TEMP, TM);

    reset(TEMP);
    a.z *= PI / 180;
    TEMP[0][0] = (float)cos(a.z);
    TEMP[0][1] = (float)sin(a.z)*(-1);
    TEMP[1][0] = (float)sin(a.z);
    TEMP[1][1] = (float)cos(a.z);
    mul(TEMP, TM);
}

void GeneralGRM(point vz, point vt) {
    reset(TEMP);
    point v1, v2;

    v1 = Cross(vt, vz);
    v2 = Cross(vz, v1);
    v1 = Normalize(v1);
    v2 = Normalize(v2);
    vz = Normalize(vz);
    TEMP[0][0] = v1.x;
    TEMP[0][1] = v1.y;
    TEMP[0][2] = v1.z;
    TEMP[1][0] = v2.x;
    TEMP[1][1] = v2.y;
    TEMP[1][2] = v2.z;
    TEMP[2][0] = vz.x;
    TEMP[2][1] = vz.y;
    TEMP[2][2] = vz.z;
}

void readFile() {
    string s = "";
    char c;
    while (c = fgetc(file)) {
        if (c == EOF)
            break;
        else {
            if (c != '\n') {
                s += c;
                continue;
            }
            else {//一個指令
                if (s.length() == 0)
                    continue;
                string in = s;
                s = "";
                if (in.length() > 0 && in[0] == '#') {//cmd line
                    cout << "\n" << in << endl;
                    continue;
                }
                else {
                    bool is_cmd = true, dot = false, neg = false;
                    string cmd = "";
                    string file_name = "";
                    int i = 0, ind = 0, dotc = 1;
                    float para[10] = { 0 }, temp = 0;
                    while (i < in.length()) {
                        if (is_cmd && in[i] != ' ') {
                            cmd += in[i];
                            i++;
                            continue;
                        }
                        if (is_cmd && in[i] == ' ') {
                            is_cmd = false;
                            i++;
                            continue;
                        }
                        if (in[i] == ' ' && in[i - 1] != ' ') {
                            temp /= dotc;//小數
                            if (neg)
                                temp *= -1;
                            para[ind] = temp;
                            ind++;
                            temp = 0;
                            dotc = 1;
                            neg = false;
                            dot = false;
                        }
                        else if (in[i] == ' ' &&in[i - 1] == ' ');//連續空白  不做事
                        else if (in[i] == '.' || in[i] == '-' || (in[i] >= '0' && in[i] <= '9')) {//number
                            if (in[i] == '.') {
                                dot = true;
                                i++;
                                continue;
                            }
                            if (in[i] == '-') {
                                neg = true;
                                i++;
                                continue;
                            }
                            temp *= 10;
                            temp += in[i] - '0';
                            if (dot)
                                dotc *= 10;
                        }
                        else {
                            while (in[i] != ' ')
                                file_name += in[i++];
                        }
                        i++;
                    }
                    temp /= dotc;//小數點
                    if (neg)
                        temp *= -1;
                    if (ind <= 9)
                        para[ind] = temp;

                    if (cmd == "")
                        continue;
                    cout << "command " << cmd << endl;
                    if (cmd == "reset")
                        reset(TM);
                    else if (cmd == "translate") {
                        point t;
                        t.x = para[0];
                        t.y = para[1];
                        t.z = para[2];
                        Translate(t, TM);
                    }
                    else if (cmd == "object") {
                        //construct data list
                        data_ptr now = new Data;
                        if (first_obj) {
                            datalist = now;
                            before = now;
                            first_obj = false;
                        }
                        before->link = now;
                        before = now;
                        now->link = 0;
                        now->O.R = para[0];
                        now->O.G = para[1];
                        now->O.B = para[2];
                        now->Kd = para[3];
                        now->Ks = para[4];
                        now->N = para[5];
                        //store object
                        if (file_name == "teapot.asc" && teapot.exist)
                            now->obj = teapot;
                        else if (file_name == "cube.asc" && cube.exist)
                            now->obj = cube;
                        else
                            now->obj = *readASC(file_name);
                        mul_obj(TM, &now->obj);
                    }
                    else if (cmd == "ambient") {
                        kaia.R = para[0];
                        kaia.G = para[1];
                        kaia.B = para[2];
                    }
                    else if (cmd == "light") {
                        light now = { 0 };
                        lightc = para[0];
                        now.lightcolor.R = para[1];
                        now.lightcolor.G = para[2];
                        now.lightcolor.B = para[3];
                        now.location.x = para[4];
                        now.location.y = para[5];
                        now.location.z = para[6];
                        lightlist[lightc - 1] = now;
                    }
                    else if (cmd == "background") {
                        bgcolor.R = para[0];
                        bgcolor.G = para[1];
                        bgcolor.B = para[2];
                    }
                    else if (cmd == "scale") {
                        point t;
                        t.x = para[0];
                        t.y = para[1];
                        t.z = para[2];
                        Scale(t);
                    }
                    else if (cmd == "rotate") {
                        point t;
                        t.x = para[0];
                        t.y = para[1];
                        t.z = para[2];
                        Rotate(t);
                    }
                    else if (cmd == "viewport") {
                        vmin.x = para[0];
                        vmin.y = para[2];
                        vmax.x = para[1];
                        vmax.y = para[3];
                        AR = (vmax.x - vmin.x) / (vmax.y - vmin.y);
                        vmax = coord2pixel(vmax, WIN_WID / 2, WIN_HEI / 2, WIN_PC);
                        vmin = coord2pixel(vmin, WIN_WID / 2, WIN_HEI / 2, WIN_PC);
                        PM[1][1] = AR;
                        pc.x = (vmax.x - vmin.x) / 2 + vmin.x;
                        pc.y = (vmax.y - vmin.y) / 2 + vmin.y;
                    }
                    else if (cmd == "observer") {
                        E.x = para[0];
                        E.y = para[1];
                        E.z = para[2];
                        COI.x = para[3];
                        COI.y = para[4];
                        COI.z = para[5];
                        tilt = para[6];
                        PMh = para[7];
                        PMy = para[8];
                        PMtheta = para[9];

                        /****EM***/
                        //T(-eye location)
                        reset(EM);
                        point t;
                        t.x = E.x*(-1);
                        t.y = E.y*(-1);
                        t.z = E.z*(-1);
                        Translate(t, EM);
                        //GRM
                        point view_vec;//v3
                        view_vec.x = COI.x - E.x;
                        view_vec.y = COI.y - E.y;
                        view_vec.z = COI.z - E.z;
                        point top_vec = { 0 };
                        top_vec.y = 1;
                        cout << "GRM" << endl;
                        GeneralGRM(view_vec, top_vec);
                        printM(TEMP);
                        mul(TEMP, EM);
                        //mirror
                        reset(TEMP);
                        TEMP[0][0] = -1;
                        cout << "Mirror" << endl;
                        printM(TEMP);
                        mul(TEMP, EM);
                        //tilt
                        reset(TEMP);
                        tilt *= (float)PI / 180;
                        TEMP[0][0] = (float)cos(tilt);
                        TEMP[0][1] = (float)sin(tilt);
                        TEMP[1][0] = (float)sin(tilt)*(-1);
                        TEMP[1][1] = (float)cos(tilt);
                        cout << "Tilt" << endl;
                        printM(TEMP);
                        mul(TEMP, EM);
                        /****PM****/
                        PMtheta *= (float)PI / 180;
                        PMtheta = tan(PMtheta);
                        reset(PM);
                        PM[1][1] = AR;
                        PM[2][2] = PMy / (PMy - PMh)*PMtheta;
                        PM[2][3] = PM[2][2] * (-1)*PMh;
                        PM[3][2] = PMtheta;
                        PM[3][3] = 0;
                    }
                    else if (cmd == "display") {
                        cout << "TM" << endl;
                        printM(TM);
                        cout << "PM" << endl;
                        printM(PM);
                        cout << "EM" << endl;
                        printM(EM);
                        glClearColor(0, 0, 0, 0);
                        glClear(GL_COLOR_BUFFER_BIT);
                        glFinish();
                        View();
                        system("pause");
                    }
                    else {//end
                        glutDestroyWindow(glutGetWindow());
                    }
                }
            }
        }
    }
    fclose(file);
}

int main(int argc, char * argv[])//这是使用glut库函数进行窗口管理
{
    system("pause");
    //set window size
    errno_t err = fopen_s(&file, argv[1], "r");
    if (err != 0)
        cout << "open file error" << endl;
    else {
        char c;
        string in = "";
        while (c = fgetc(file)) {
            if (c == '\n')
                break;
            else
                in += c;
        }
        int temp = 0;
        for (int i = 0; i < in.length(); i++) {
            if (in[i] == ' ') {
                WIN_WID = temp;
                temp = 0;
                continue;
            }
            temp *= 10;
            temp += in[i] - '0';
        }
        WIN_HEI = temp;
        WIN_PC.x = floor(WIN_WID / 2 + 0.5);
        WIN_PC.y = floor(WIN_HEI / 2 + 0.5);
    }

    glutInit(&argc, argv);//使用glut库需要进行初始化
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);//设定窗口显示模式，颜色模型和缓存，这里是RGB颜色模型和单缓存
    glutInitWindowPosition(10, 10);//设定窗口的初始位置，屏幕左上角为原点，单位为像素
    glutInitWindowSize(WIN_WID, WIN_HEI);//设定窗口的大小
    glutCreateWindow("HW4");//创建一个窗口，参数是窗口标题名
    glutDisplayFunc(myDisplay);//将myDisplay指定为当前窗口的显示内容函数
    Initial();
    readFile();
    glutMainLoop();//使窗口框架运行起来，使显示回调函数开始工作
    return 0;
}