//Liu Elina Graham Scan Hull

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <chrono>
#include <stack>
#include <iomanip> //added to use setprecision() for output vals
#include <algorithm> //to use sort()
const int DIM = 400;
int ARR[DIM][DIM];
using namespace std;

class Point{
    private:
        double cx = 0.0;
        double cy = 0.0;
        int x = scale(cx);
        int y = scale(cy);

        static int scale(double val){
            int scaled = val * DIM;
            return scaled;
        }
    public:
        //form constructors
        Point(){
            cx = random_val();
            cy = random_val();
        }
        Point(double one, double two){
            cx = one;
            cy = two;
        }
        double random_val(){ //generates a random
            double val = (double)rand()/(double)RAND_MAX;
            return val;
        }
        //make access to private vars
        double getcx(){
            return cx;
        }
        double getcy(){
            return cy;
        }
        void setcx(double xval){
            cx = xval;
        }
        void setcy(double yval){
            cy = yval;
        }
        double distoLine(Point p1, Point p2){
            return ((p2.getcx() - p1.getcx())*(p1.getcy()-cy)- (p1.getcx()-cx)*(p2.getcy()-p1.getcy())) / (sqrt(pow(p1.getcx()-p2.getcx(),2) + pow(p1.getcy()-p2.getcy(),2)));
        }
};

///// rendering into ppm
bool within_board(int x, int y){
    if(0>x || x>DIM){
        return false;
    }
    if(0>y || y>DIM){
        return false;
    }
    else{
        return true;
    }
}
void drawCircle(int xc, int yc, int x, int y)
{
    if(within_board(yc+y, xc+x)){
        ARR[yc+y][xc+x] = 0;
    }
    if(within_board(yc-y, xc+x)){
        ARR[yc-y][xc+x] = 0;
    }
    if(within_board(yc-y, xc-x)){
        ARR[yc-y][xc-x] = 0;
    }
    if(within_board(yc+y, xc-x)){
        ARR[yc+y][xc-x] = 0;
    }
    if(within_board(yc+x, xc-y)){
        ARR[yc+x][xc-y] = 0;
    }
    if(within_board(yc+x, xc+y)){
        ARR[yc+x][xc+y] = 0; 
    }
    if(within_board(yc-x, xc-y)){
        ARR[yc-x][xc-y] = 0;
    }
    if(within_board(yc-x,xc+y)){
        ARR[yc-x][xc+y] = 0;
    }
}
void draw_circle(int xc, int yc, int r){ //h, k, radius
    int x = 0, y = r;
    int d = 3 - 2 * r;
    drawCircle(xc, yc, x, y);
    while (y >= x)
    {
        x++;
        if (d > 0)
        {
            y--;
            d = d + 4 * (x - y) + 10;
        }
        else{
            d = d + 4 * x + 6;
        }
        drawCircle(xc, yc, x, y);
    }
}
void bresenham(int x1, int y1, int x2, int y2)
{
    if(abs(double(y2*1.0-y1)/(x2-x1)) > 1){ //da is y
    //for DA=y-axis
        if(y1>y2){ //swap the points if right point given before left
            int temp = x2;
            x2=x1;
            x1=temp;
            temp =y2;
            y2=y1;
            y1=temp;
        }
        int deltax = abs(x2-x1);
        int deltay = y2-y1;
        int j = x1;
        int e = deltax-deltay; //always negative
        bool xdir = (x2>x1); //true means y ascneding
        for(int i = y1; i<=y2; i++){
            //draw the pixel (i,j)
            if (within_board(i,j)){
                ARR[i][j] = 0; //color in the 2d array point
            }
            if(e>=0){
                if(xdir){
                    j++;
                }
                else{
                    j--;
                }
                e-=deltay;
            }
            e+=deltax;
        }
    }
    else{
     //for DA=x-axis (slope domain: -1,1)
        if(x1>x2){ //swap the points if right point given before left
            int temp = x2;
            x2=x1;
            x1=temp;
            temp =y2;
            y2=y1;
            y1=temp;
        }
        int deltax = x2-x1;
        int deltay = abs(y2-y1);
        int j = y1;
        int e = deltay-deltax; //always negative
        bool ydir = (y2>y1); //true means y ascneding
        for(int i = x1; i<=x2; i++){
            //draw the pixel (i,j)
            if (within_board(j,i)){
                ARR[j][i] = 0; //color in the 2d array point
            }
            if(e>=0){
                if(ydir){
                    j++;
                }
                else{
                    j--;
                }
                e-=deltax;
            }
            e+=deltay;
        }
    }
}
//arraythings
void ppm_it(ofstream &f){//prints the ppm
    f<<"P3 "<< DIM<< " "<<DIM<< " 1\n";
    for(int i=DIM-1;i>=0;i--){
        for(int j=0;j<DIM;j++){
            int val = ARR[i][j];
            f<<val<<" "<<val<<" "<<val<<" ";
        }
        f<<"\n";
    }
}

void fill_array(){
    for(int i =0;i<DIM;i++){
        for(int j=0;j<DIM;j++){
            ARR[i][j]=1;
        }
    }
}
int scale(double val){
    int scaled = val * DIM;
    return scaled;
}
/////


void part0(int numpts){ //generates numpts number of points
    srand ( time(0) );
    ofstream outfile;
    outfile.open("points.txt");
    for(int i = 0; i<numpts;i++){ //generate numpts random points
        Point r;
        outfile<<fixed<<setprecision(23)<<r.getcx()<<"  "<<r.getcy()<<"\n";
    }
    outfile.close();
}
void findHull(vector<Point> &cvHull, vector<Point> &Sk, Point &P, Point &Q){
    if(Sk.size()==0){
        return;
    }
    double maxdis = 0;
    Point C;
    for(unsigned int i = 0; i<Sk.size();i++){ //find pt C farthest from PQ on the right side
        double distance = Sk[i].distoLine(P, Q);
        if(abs(distance) > maxdis){
            maxdis = abs(distance);
            C = Sk[i];
        }
    }
    //find index of P
    int p_ind; 
    for(unsigned int i = 0; i<cvHull.size(); i++){
        if(cvHull[i].getcx()==P.getcx()  && cvHull[i].getcy()==P.getcy()){ //if those pts are equal
            p_ind = i;
            break;
        }
    }
    cvHull.insert(cvHull.begin() + p_ind + 1, C);
    
    //make s1 and s2
    vector<Point> rightPC;
    vector<Point> rightCQ;
    
    for(unsigned int i = 0; i<Sk.size(); i++){
        //sort by right and left side of line
        if(Sk[i].distoLine(P, C)>0){ //right of PC
            rightPC.push_back(Sk[i]);
        }
        else if (Sk[i].distoLine(C, Q)>0){
            rightCQ.push_back(Sk[i]);
        }    
    }
    
    findHull(cvHull, rightPC, P, C);
    findHull(cvHull, rightCQ, C, Q);
    return;
}

vector<Point> QuickHull(vector<Point> &points){ //n is number of points
    //find convex hull from set points 
    //find left and right most pts 
    vector<Point> cvHull;
    sort(points.begin(), points.end(), [](Point lhs, Point rhs ) 
    {
        return lhs.getcx() < rhs.getcx();
    });
    
    Point A = points[0];
    Point B = points[points.size()-1];
    cvHull.push_back(A);
    cvHull.push_back(B);
    
    vector<Point> rightsAB;
    vector<Point> rightsBA;
    
    for(unsigned int i = 0; i<points.size(); i++){
        //sort by right and left side of line
        if(points[i].distoLine(A, B)<=0){ //left of AB
            rightsBA.push_back(points[i]);
        }
        else{
            rightsAB.push_back(points[i]);
        }    
    }
    findHull(cvHull, rightsAB, A, B);
    findHull(cvHull, rightsBA, B, A);
    cvHull.push_back(A);
    return cvHull;
    
}

void printHull1(vector<Point> &hull, vector<Point> &points, ofstream &f){
  //  cout<<"Points in Convex Hull: \n";
    fill_array();
//     for(unsigned int i  = 0; i<hull.size();i++){
//         cout<<"("<<hull[i].getcx()<<", "<<hull[i].getcy()<<") , ";
//     }
//    cout<<"\nAll Points: \n";
    for(unsigned int i = 0; i<points.size(); i++){
        draw_circle(scale(points[i].getcx()), scale(points[i].getcy()), 2); //draw each pt
       // cout<<"("<<points[i].getcx()<<", "<<points[i].getcy()<<") , ";
    }
    
    //hull
    for(unsigned int i = 0; i<hull.size()-1;i++){ //connects lines of hull
        Point a = hull[i];
        Point b = hull[i+1];
        bresenham(scale(a.getcx()), scale(a.getcy()), scale(b.getcx()), scale(b.getcy()));
    }
    ppm_it(f);
    
}
void part1(){
    part0(60); //generate 60 pts, save into points
    vector<Point> points;
    ifstream infile;
    infile.open("points.txt");
    ofstream outppm;
    outppm.open("quickhull.ppm");
    double x,y;
    while(infile >>x>>y){ //read in all points
        Point a(x, y); 
        points.push_back(a);
    }
    infile.close();
    vector<Point> convexHull;
    convexHull = QuickHull(points);
    printHull1(convexHull, points, outppm);
    outppm.close();
}


/////part2
int turntype(Point a, Point b, Point c){ //returns orientation of turn from AB to C
    double val = (b.getcy()- a.getcy()) * (c.getcx() -b.getcx()) - (b.getcx() -a.getcx()) * (c.getcy()- b.getcy());
    if (val == 0.0){// collinear
        return 0;
    } 
    else if (val>0){ //clockwise, right turn
        return 1;
    }
    else{ //counterclockwise, left turn
        return 2;
    }
}
double dotcos(Point a, Point b){ //finds costheta, by dotting a and b
    Point vect1(1-a.getcx(),0); //vector a to 1,a.cy
    Point vect2(b.getcx()-a.getcx(),b.getcy()-a.getcy()); //vector ab
    
    double dot = vect1.getcx()*vect2.getcx()+vect1.getcy()*vect2.getcy();
    double magn = sqrt(vect1.getcx()*vect1.getcx()+vect1.getcy()*vect1.getcy()) * sqrt(pow(vect2.getcx(),2)+pow(vect2.getcy(),2));
    return dot/magn;
}
double distSq(Point p1, Point p2){
    return (p1.getcx() - p2.getcx())*(p1.getcx() - p2.getcx()) + (p1.getcy() - p2.getcy())*(p1.getcy() - p2.getcy());
}
stack<Point> GrahamScan(vector<Point> &points, stack<Point> &cvHull){
    sort(points.begin(), points.end(), [](Point lhs, Point rhs ) //find lowest y-coord
    {
        return lhs.getcy() < rhs.getcy();
    });
    
    Point p0 = points[0];
    
    sort(points.begin()+1, points.end(), [p0](Point lhs, Point rhs ) //sort by polar angle using cos to approx
    {
        //return true if lhs should be before rhs
        //cos decreasing
        double coslhs = dotcos(p0, lhs);
        double cosrhs = dotcos(p0, rhs);
        if (cosrhs==coslhs){//same polar angle
            return distSq(p0, lhs) > distSq(p0, rhs); //consider farthest from p0
        } 
        return coslhs > cosrhs;
    }); 
    
    cvHull.push(p0);
    cvHull.push(points[1]);
    
    for(unsigned int i = 2; i<points.size();i++){
        Point b = cvHull.top();
        cvHull.pop();
        Point a = cvHull.top();
        while(cvHull.size()>1 && turntype(a, b, points[i]) < 2){ //if the turn is clockwise, remove points
            b = cvHull.top();
            cvHull.pop();
            a = cvHull.top();
        }
        cvHull.push(b);
        cvHull.push(points[i]);
    }
    cvHull.push(p0); //add first point again to connect the hull
    return cvHull;
}
void printHull2(vector<Point> &points, stack<Point> &hull, ofstream &f){
    fill_array();
    stack<Point> copyHull = hull;
   // cout<<"Points in Hull: \n";
//     cout<<copyHull.size()-1<<"\n";
    while(copyHull.size()>1){
        Point b = copyHull.top();
        copyHull.pop();
        Point a = copyHull.top();
//         cout<<"("<<b.getcx()<<", "<<b.getcy()<<") , ";
        bresenham(scale(a.getcx()), scale(a.getcy()), scale(b.getcx()), scale(b.getcy()));
    }
//     cout<<"\nAll Points: \n";
    for(unsigned int i = 0; i<points.size(); i++){
        draw_circle(scale(points[i].getcx()), scale(points[i].getcy()), 3); //draw each pt
//         cout<<"("<<points[i].getcx()<<", "<<points[i].getcy()<<") , ";
    }
//     cout<<"\n";
    ppm_it(f);
}

void part2(){
    part0(60);
    vector<Point> points;
    ifstream infile;
    infile.open("points.txt");
    ofstream outppm;
    outppm.open("grahamscan.ppm");
    double x,y;
    while(infile >>x>>y){ //read in all points
        Point a(x, y); 
        points.push_back(a);
    }
    infile.close();
    stack<Point> convexHull;
    convexHull = GrahamScan(points, convexHull);
    printHull2(points, convexHull, outppm);
    outppm.close();
}

int main(){
    part1();
    part2();
}