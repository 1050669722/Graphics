//
// Created by ASUS on 2023/5/22.
//

#include <iostream>

#include <bits/stdc++.h>
using namespace std;

struct Point
{
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};

typedef Point Vector;

Vector operator+(Vector A, Vector B)
{
    return Vector(A.x + B.x, A.y + B.y);
}

Vector operator-(Point A, Point B)
{
    return Vector(A.x - B.x, A.y - B.y);
}

Vector operator*(Vector A, double p)
{
    return Vector(A.x * p, A.y * p);
}

Vector operator/(Vector A, double p)
{
    return Vector(A.x / p, A.y / p);
}

const double eps = 1e-10;
int dcmp(double x)
{
    if (fabs(x) < eps) return 0;
    else return x < 0 ? -1 : 1;
}

bool operator<(const Point& a, const Point& b)
{
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

bool operator==(const Point& a, const Point& b)
{
    return dcmp(a.x - b.x) == 0 && dcmp(a.y - b.y) == 0;
}

double Dot(Vector A, Vector B)
{
    return A.x * B.x + A.y * B.y;
}

double Cross(Vector A, Vector B)
{
    return A.x * B.y - A.y * B.x;
}

double Length(Vector A)
{
    return sqrt(Dot(A, A));
}

double Angle(Vector A, Vector B)
{
    return acos(Dot(A, B) / Length(A) / Length(B));
}

Vector Rotate(Vector A, double rad)
{
    return Vector(A.x * cos(rad) - A.y * sin(rad), A.x * sin(rad) + A.y * cos(rad));
}

Vector Normal(Vector A)
{
    double L = Length(A);
    return Vector(-A.y / L, A.x / L);
}

Point GetLineIntersection(Point P, Vector v, Point Q, Vector w)
{
    Vector u = P - Q;
    double t = Cross(w, u) / Cross(v, w);
    return P + v * t;
}

double DistanceToLine(Point P, Point A, Point B)
{
    Vector v1 = B - A, v2 = P - A;
    return fabs(Cross(v1, v2)) / Length(v1);
}

double DistanceToSegment(Point P, Point A, Point B)
{
    if (A == B) return Length(P - A);
    Vector v1 = B - A, v2 = P - A, v3 = P - B;
    if (dcmp(Dot(v1, v2)) < 0) return Length(v2);
    else if (dcmp(Dot(v1, v3)) > 0) return Length(v3);
    else return fabs(Cross(v1, v2)) / Length(v1);
}

Point GetLineProjection(Point P, Point A, Point B)
{
    Vector v = B - A;
    return A + v * (Dot(v, P - A) / Dot(v, v));
}

bool SegmentProperIntersection(Point a1, Point a2, Point b1, Point b2)
{
    double c1 = Cross(a2 - a1, b1 - a1), c2 = Cross(a2 - a1, b2 - a1),
            c3 = Cross(b2 - b1, a1 - b1), c4 = Cross(b2 - b1, a2 - b1);
    return dcmp(c1) * dcmp(c2) < 0 && dcmp(c3) * dcmp(c4) < 0;
}

bool OnSegment(Point p, Point a1, Point a2)
{
    return dcmp(Cross(a1 - p, a2 - p)) == 0 && dcmp(Dot(a1 - p, a2 - p)) < 0;
}

int ConvexHull(Point* p, int n, Point* ch)
{
    sort(p, p + n);
    int m = 0;
    for (int i = 0; i < n; i++) {
        while (m > 1 && Cross(ch[m - 1] - ch[m - 2], p[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = p[i];
    }
    int k = m;
    for (int i = n - 2; i >= 0; i--) {
        while (m > k && Cross(ch[m - 1] - ch[m - 2], p[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = p[i];
    }
    if (n > 1) m--;
    return m;
}

int main(int argc, char* argv[])
{
    std::cout << "Hello World" << std::endl;

    // int n;
    // scanf("%d", &n);
    // Point p[n], ch[n];
    // for (int i = 0; i < n; i++) {
    //     scanf("%lf%lf", &p[i].x, &p[i].y);
    // }

    int n = 6;
    Point ch[n];
    Point p[]
    {
        Point(0, 0), Point(1, 2), Point(2, 0),
        Point(1, 1), Point(0.5, 1), Point(1.5, 1),
    };
    // Point p[n];
    // p[0] = Point(0, 0);
    // p[1] = Point(1, 2);
    // p[2] = Point(2, 0);
    // p[3] = Point(1, 1);
    // p[4] = Point(0.5, 1);
    // p[5] = Point(1.5, 1);
    // for (int i = 0; i < n; ++i)
    // {
    //     std::cout << p[i].x << ", " << p[i].y << std::endl;
    // }

    int m = ConvexHull(p, n, ch);
    printf("Convex Hull:\n");
    for (int i = 0; i < m; i++) {
        printf("(%.1lf,%.1lf)\n", ch[i].x, ch[i].y);
    }

    return 0;
}