#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <iomanip>
#include <algorithm>
#include <fstream>

double max(double a, double b)
{
    return a > b ? a : b;
}

class Vector3
{
public:
    double x, y, z;

    Vector3() {}

    Vector3(const Vector3 &a) : x(a.x), y(a.y), z(a.z) {}
    Vector3(double nx, double ny, double nz) : x(nx), y(ny), z(nz) {}

    Vector3 operator+(Vector3 const &vec)
    {
        Vector3 res;
        res.x = x + vec.x;
        res.y = y + vec.y;
        res.z = z + vec.z;
        return res;
    }

    Vector3 operator-(Vector3 const &vec)
    {
        Vector3 res;
        res.x = x - vec.x;
        res.y = y - vec.y;
        res.z = z - vec.z;
        return res;
    }

    double dot(Vector3 const &vec)
    {
        return x * vec.x + y * vec.y + z * vec.z;
    }

    Vector3 times(double a)
    {
        Vector3 res;
        res.x = x * a;
        res.y = y * a;
        res.z = z * a;
        return res;
    }

    double dist(Vector3 const &vec)
    {
        return sqrt(pow(x - vec.x, 2) + pow(y - vec.y, 2) + pow(z - vec.z, 2));
    }

    void print()
    {
        std::cout << "x: " << x << " "
                  << "y: " << y << " "
                  << "z: " << z << " "
                  << std::endl;
    }
};

class Particle
{
public:
    Vector3 position;
    Vector3 velocity;
    Vector3 acc;
    Vector3 prevAcc;
    Vector3 force;
    double density;
    double pressure;
};

double h = 0.22;                                                       // 影響範囲
double particleMass = 0.0001;                                          // 粒子の質量
double densityCoef = particleMass * 315.0 / (64.0 * M_PI * pow(h, 9)); // 密度計算で用いる

// 密度の計算
void CalcDensity(std::vector<Particle> &particles)
{
    double h2 = pow(h, 2);
    for (int i = 0; i < particles.size(); i++)
    {
        Particle *nowParticle = &particles[i];
        double sum = 0;
        for (int j = 0; j < particles.size(); j++)
        {
            if (i == j)
            {
                continue;
            }
            Particle nearParticle = particles[j];
            Vector3 diff = nearParticle.position - nowParticle->position;
            double r2 = diff.dot(diff);

            if (r2 < h2)
            {

                double c = h2 - r2;
                sum += pow(c, 3);
            }
        }
        // 密度
        nowParticle->density = sum * densityCoef + 0.000001;
    }
}

double pressureStiffness = 10.0; // 圧力係数k
double restDensity = 0;          // 静止密度ρ0(ローゼロ)

// 圧力
void calcPressure(std::vector<Particle> &particles)
{
    for (int i = 0; i < particles.size(); i++)
    {
        particles[i].pressure = pressureStiffness * (particles[i].density - restDensity);
    }
}

double pressureCoef = particleMass * 45.0 / (M_PI * pow(h, 6)); // 圧力項の計算で用いる

// 圧力項の計算
Vector3 calcPressurePart(std::vector<Particle> &particles, int target)
{
    double h2 = pow(h, 2);
    Particle nowParticle = particles[target];
    Vector3 gradSum = Vector3(0.0, 0.0, 0.0);
    for (int neighbor = 0; neighbor < particles.size(); neighbor++)
    {
        if (target == neighbor)
        {
            continue;
        }

        Particle nearParticle = particles[neighbor];
        Vector3 diff = nearParticle.position - nowParticle.position;
        double r2 = diff.dot(diff);

        if (r2 < h2)
        {
            double scalar = pow(h - sqrt(r2), 2.0) * pressureCoef * (nearParticle.pressure - nowParticle.pressure) / (2.0 * nearParticle.density);
            Vector3 vectorNormalized = diff.times(1.0 / sqrt(r2));
            gradSum = gradSum + vectorNormalized.times(scalar);
        }
    }
    return gradSum.times(-1.0 / nowParticle.density);
}

double viscosityStiffness = 0.8;                                                      // 粘性係数v
double viscosityCoef = viscosityStiffness * particleMass * 45.0 / (M_PI * pow(h, 6)); // 粘性項の計算で用いる

// 粘性項の計算
Vector3
calcViscosityPart(std::vector<Particle> &particles, int target)
{
    double h2 = pow(h, 2);
    Particle nowParticle = particles[target];
    Vector3 sum = Vector3(0.0, 0.0, 0.0);
    for (int neighbor = 0; neighbor < particles.size(); neighbor++)
    {
        if (target == neighbor)
        {
            continue;
        }
        Particle nearParticle = particles[neighbor];
        Vector3 diff = nearParticle.position - nowParticle.position;
        double r2 = diff.dot(diff);

        if (r2 < h2)
        {
            double scalar = viscosityCoef * (h - sqrt(r2)) / nearParticle.density;

            sum = sum + (nearParticle.velocity - nowParticle.velocity).times(scalar);
        }
    }
    return sum;
}

// 各壁の法線ベクトル
Vector3 wallNorthN = Vector3(0, -1, 0);
Vector3 wallSouthN = Vector3(0, 1, 0);
Vector3 wallEastN = Vector3(-1, 0, 0);
Vector3 wallWestN = Vector3(1, 0, 0);
Vector3 wallBottomN = Vector3(0, 0, 1);
Vector3 wallTopN = Vector3(0, 0, -1);

// 各壁を表す平面 wallNorthの場合はy=0.5の平面ということ
Vector3 wallNorth = Vector3(0, 0.5, 0);
Vector3 wallSouth = Vector3(0, -0.5, 0);
Vector3 wallEast = Vector3(0.5, 0, 0);
Vector3 wallWest = Vector3(-0.5, 0, 0);
Vector3 wallBottom = Vector3(0, 0, 0);
Vector3 wallTop = Vector3(0, 0, 1);

double spring = 10.0; // バネ係数
double damper = 20;   // ダンパ係数

// 衝突判定(ペナルティ法)
// 壁から跳ね返す力の合力を求める
Vector3
calcCollisionForce(std::vector<Particle> &particles, int target)
{
    Particle nowParticle = particles[target];
    Vector3 f = Vector3(0, 0, 0);
    //  北側の壁からの粒子を跳ね返す力
    if (nowParticle.position.y - wallNorth.y > 0)
    {
        double penetrate = nowParticle.position.y - wallNorth.y;
        Vector3 fNorth = wallNorthN.times(spring * penetrate - damper * nowParticle.velocity.dot(wallNorthN));
        f = f + fNorth;
    }

    // // 南側の壁からの粒子を跳ね返す力
    if (nowParticle.position.y - wallSouth.y < 0)
    {
        double penetrate = wallSouth.y - nowParticle.position.y;
        Vector3 fSouth = wallSouthN.times(spring * penetrate - damper * nowParticle.velocity.dot(wallSouthN));
        f = f + fSouth;
    }

    // // 東側の壁からの粒子を跳ね返す力
    if (nowParticle.position.x - wallEast.x > 0)
    {
        double penetrate = nowParticle.position.x - wallEast.x;
        Vector3 fEast = wallEastN.times(spring * penetrate - damper * nowParticle.velocity.dot(wallEastN));
        f = f + fEast;
    }

    // // 西側の壁からの粒子を跳ね返す力
    if (nowParticle.position.x - wallWest.x < 0)
    {
        double penetrate = wallWest.x - nowParticle.position.x;
        Vector3 fWest = wallWestN.times(spring * penetrate - damper * nowParticle.velocity.dot(wallWestN));
        f = f + fWest;
    }

    // 下側の壁からの粒子を跳ね返す力
    if (nowParticle.position.z - wallBottom.z < 0)
    {
        double penetrate = wallBottomN.z - nowParticle.position.z;
        Vector3 fBottom = wallBottomN.times(spring * penetrate - damper * nowParticle.velocity.dot(wallBottomN));
        f = f + fBottom;
    }

    // 上側の壁からの粒子を跳ね返す力

    return f;
}

void simulate(std::vector<Particle> &particles, double deltaTime)
{
    std::ofstream file;
    file.open("fluid.csv");

    for (int count = 0; count < 1000; count++)
    {
        // 密度・圧力の値を更新
        CalcDensity(particles);
        calcPressure(particles);

        // std::cout << "pressure: " << particles[0].pressure << std::endl;
        // std::cout << "density: " << particles[0].density << std::endl;
        // std::cout << "pos z: " << particles[0].position.z << std::endl;
        std::cout << "pos acc" << std::endl;
        particles[0].position.print();
        particles[0].acc.print();

        for (int targetParticle = 0; targetParticle < particles.size(); targetParticle++)
        {

            Particle *nowParticle = &particles[targetParticle];

            // 圧力項 + 粘性項 + 衝突による力 + 重力
            nowParticle->acc = calcPressurePart(particles, targetParticle) + calcViscosityPart(particles, targetParticle) + calcCollisionForce(particles, targetParticle) + Vector3(0, 0, -0.6);

            // リープ・フロッグ法によるvelocity, positionの計算
            nowParticle->velocity = nowParticle->velocity + (nowParticle->prevAcc + nowParticle->acc).times(0.5 * deltaTime);

            nowParticle->position = nowParticle->position + nowParticle->velocity.times(deltaTime) + nowParticle->acc.times(0.5 * pow(deltaTime, 2));
            file << targetParticle << "," << nowParticle->position.x << "," << nowParticle->position.y << "," << nowParticle->position.z << "\n";

            nowParticle->prevAcc = nowParticle->acc;
        }
    }

    file.close();
}

int main()
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(-0.01, 0.05);

    int pCount = 20000;
    std::vector<Particle> particles(pCount);
    for (int i = 0; i < particles.size(); i++)
    {

        particles[i].velocity = Vector3(0.0, 0.0, 0.0);
        particles[i].position = Vector3(distr(eng), distr(eng), 0.9 + distr(eng));
        particles[i].acc = Vector3(0.0, 0.0, 0.0);
        particles[i].prevAcc = Vector3(0.0, 0.0, 0.0);
    }

    simulate(particles, 0.05);
    return 0;
}