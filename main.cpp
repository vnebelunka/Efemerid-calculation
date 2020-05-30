
#include <iostream>
#include <cstdio>
#include <cmath>
#define PI 3.141592653589793238
const double eps=0.00000001;
const double ekl= 23.44760833/180*PI;//наклон земного экватора к еклиптике
using namespace std;
///рассчет эпохи по строке
int TimeToUlian(int year,int month,int day)
    {
       int a=(14-month)/12,y=year+4800-a,m=month+12*a-3;
       int UlianTime=day+(153*m+2)/5+365*y+y/4-y/100+y/400-32045;
       return UlianTime;
    }

///функция решает уравнение Кеплера
double Find_E(double M,double e)
{
    double E0,E1=0;
    E0=M;
    while(fabs(E0-E1)>eps)
    {
        E1=E0;
        E0=e*sin(E1)+M;
    }
    return E1;
}

///время в 3 форматах: день месяц год, сокращенная эпоха, юлианский день
struct TIME
{
    char Epoh[5];
    int day,month, year;
    long long UlianTime;
    void EpohToTime()
    {
    int preyear;
    switch(Epoh[0])
    {
       case 'I':
        preyear=18;
        break;
    case 'J':
        preyear=19;
        break;
    case 'K':
        preyear=20;
        break;
    }
    year=preyear*100+(Epoh[1]-'0')*10+Epoh[2]-'0';
    if(Epoh[3]>='0'&&Epoh[3]<='9') month=Epoh[3]-'0';
    else month=Epoh[3]-'A'+10;
    if(Epoh[4]>='0'&&Epoh[4]<='9') day=Epoh[4]-'0';
    else day=Epoh[4]-'A'+10;
    }
};

///характеристика орбиты планеты
struct PLANET
{
        double e, //экцентриситет планеты
           M0,//аномалия на эпоху планеты
           M,//средняя аномалия в момент t
           i, //наклонение планеты
           w,//аргумент перицентра планеты
           W,//долгота восходящего узла планеты
           n,//среднее движение планеты
           a;//большая полуось
        TIME epoh;//эпоха
};
struct COORD
{
    double x,y,z;
};

///расчет характеристик обриты земли на момент t
COORD Earth_input(long long t,PLANET Earth)
{
    Earth.a=1.00000101778;
    Earth.e= 0.0167086342;
    double T=(t-2451545.0)/365250.0;
    Earth.i=469.97289/648000*PI*T;
   long double delta=100.46645683/180*PI+1295977422.83429/648000*PI*T;
   long double per=102.93734808/180*PI+ 11612.35290/648000*PI*T;
    Earth.W= 174.87317577/180*PI- 8679.27034/648000*PI*T;
    Earth.M=delta-per;
    Earth.w=per-Earth.W;
    double E=Find_E(Earth.M,Earth.e),
    check4=E-Earth.e*sin(E)-Earth.M,//0
    //орбитальная система координат
    osc1=Earth.a*(cos(E)-Earth.e),
    osc2=Earth.a*sqrt(1-Earth.e*Earth.e)*sin(E),
    //коэффы перехода от  орбитальной к гео-эклиптической
    Px=cos(Earth.W)*cos(Earth.w)-sin(Earth.w)*sin(Earth.W)*cos(Earth.i),
    Py=cos(Earth.w)*sin(Earth.W)+sin(Earth.w)*cos(Earth.W)*cos(Earth.i),
    Pz=sin(Earth.w)*sin(Earth.i),
    Qx=-sin(Earth.w)*cos(Earth.W)-cos(Earth.w)*sin(Earth.W)*cos(Earth.i),
    Qy=-sin(Earth.w)*sin(Earth.W)+cos(Earth.w)*cos(Earth.W)*cos(Earth.i),
    Qz=cos(Earth.w)*sin(Earth.i),
    check1=Px*Px+Py*Py+Pz*Pz,//1
    check2=Qx*Qx+Qz*Qz*Qy*Qy,//1
    check3=Px*Qx+Pz*Qz+Py*Qy,//0
    //гелиоцентрические эклиптические координаты Земли
    x=Px*osc1+Qx*osc2,git push -u origin master
    y=Py*osc1+Qy*osc2,
    z=Pz*osc1+Qz*osc2;
    COORD EarthCoordekl,EarthCoordekvat;
    //гелиоцентрические эклиптические координаты Земли
    EarthCoordekl.x=x;EarthCoordekl.y=y;EarthCoordekl.z=z;
    //гелиоцентрические экваториальные координаты Земли
    EarthCoordekvat.x=EarthCoordekl.x;
    EarthCoordekvat.y=EarthCoordekl.y*cos(ekl)-EarthCoordekl.z*sin(ekl);
    EarthCoordekvat.z=EarthCoordekl.y*sin(ekl)-EarthCoordekl.z*cos(ekl);
    return EarthCoordekvat;
}
///гелиоцентрические эклиптические координаты планеты на время t
COORD FindCoordSunEklipt(PLANET planet,long long t)
{
    planet.M=planet.M0+planet.n*(t-planet.epoh.UlianTime);
    double E=Find_E(planet.M,planet.e),
    osc1=planet.a*(cos(E)-planet.e),//орбитальная система координат
    osc2=planet.a*sqrt(1-planet.e*planet.e)*sin(E),//орбитальная система координат
    //коэффы перехода от  орбитальной к гео-эклиптической
    Px=cos(planet.W)*cos(planet.w)-sin(planet.w)*sin(planet.W)*cos(planet.i),
    Py=cos(planet.w)*sin(planet.W)+sin(planet.w)*cos(planet.W)*cos(planet.i),
    Pz=sin(planet.w)*sin(planet.i),
    Qx=-sin(planet.w)*cos(planet.W)-cos(planet.w)*sin(planet.W)*cos(planet.i),
    Qy=-sin(planet.w)*sin(planet.W)+cos(planet.w)*cos(planet.W)*cos(planet.i),
    Qz=cos(planet.w)*sin(planet.i);
    COORD coordSun;
    //гелиоцентрические эклиптические координаты планеты
    coordSun.x=Px*osc1+Qx*osc2;
    coordSun.y=Py*osc1+Qy*osc2;
    coordSun.z=Pz*osc1+Qz*osc2;
    return coordSun;
}
///расчет гелиоцентрических экваториальных координат планеты
///по гелиоцентрическим эклиптическим
COORD FindCoordSunEkvat(COORD CoordSun)
{
    COORD Coord;
    Coord.x=CoordSun.x;
    Coord.y=CoordSun.y*cos(ekl)-CoordSun.z*sin(ekl);
    Coord.z=CoordSun.y*sin(ekl)+CoordSun.z*cos(ekl);
    return Coord;
}
///геоцентрические координаты планеты по гелицентрическим экваториальным
///и геоцентрическим координатам солнца
COORD FindCoordEarth(COORD CoordSunEkvat,COORD GeoSUN)
{
    COORD coord;
    coord.x=CoordSunEkvat.x+GeoSUN.x;
     coord.y=CoordSunEkvat.y+GeoSUN.y;
      coord.z=CoordSunEkvat.z+GeoSUN.z;
      return coord;
}
///ввод информации об орбите искомой планеты
void InputPlanet(FILE *fd,PLANET &planet)
{
    fscanf(fd,"Epoh = %5s\nM = %lf\nPeri(w) = %lf\nNode(W) = %lf\nIncl(I) = %lf\ne = %lf\nn = %lf\na = %lf",
           &planet.epoh.Epoh,
           &planet.M0,
           &planet.w,
           &planet.W,
           &planet.i,
           &planet.e,
           &planet.n,
           &planet.a);
    planet.epoh.EpohToTime();
    TimeToUlian(planet.epoh.year,planet.epoh.month,planet.epoh.day);
     planet.epoh.EpohToTime();
    planet.epoh.UlianTime=TimeToUlian(planet.epoh.year,planet.epoh.month,planet.epoh.day);
    planet.i=planet.i*PI/180;
    planet.w=planet.w*PI/180;
    planet.W=planet.W*PI/180;
    planet.M0=planet.M0*PI/180;
    planet.n=planet.n*PI/180;
}
void AlphaClock(double AlphaGrad,FILE *fd1)
{
    double AlphaSek=AlphaGrad*240;
    int AlphaHour=floor(AlphaSek/(3600));
    AlphaSek=AlphaSek-AlphaHour*3600;
    int AlphaMin=floor(AlphaSek/60);
    AlphaSek-=AlphaMin*60;
    fprintf(fd1,"Alpha = %dh%dm%.2lfs\n",AlphaHour,AlphaMin,AlphaSek);
}
void SigmaGrad(double sigma,FILE *fd1)
{
    double SigmaSek=sigma*3600;
    int SigmaGrad=trunc(SigmaSek/3600);
    SigmaSek-=SigmaGrad*3600;
    int SigmaMin=trunc(SigmaSek/60);
    SigmaSek-=SigmaMin*60;
    if(SigmaMin<0)SigmaMin*=-1;
    if(SigmaSek<0)SigmaSek*=-1;
    fprintf(fd1,"Sigma = %d°%d'%.2lf\"\n",SigmaGrad,SigmaMin,SigmaSek);
}
int main()
{
    FILE *fd=fopen("input.txt","r");
    FILE *fd1=fopen("output.txt","w");
    PLANET planet,Earth;
    int tMonth,tDay,tYear,//время, на которое расчитывается эфемерида
                    ulianT;//юлианский день входного времени
    COORD GeoSun,//геоцентрические координаты солнца на t
          GelioEarth,//гелиоцентрические координаты Земли на t
          CoordPlanetSunEklipt,//гелиоцентрические эклиптечиские координаты планеты
          CoordPlanetSunEkvat,//гелиоцентрические экваториальные координаты планеты
          CoordEarth;//геоцентрические координаты планеты
    double alpha,sigma;//прямое восхождение и склонение.
    InputPlanet(fd,planet);
    fscanf(fd,"\nDate of efemirid(yyyy mm dd):\n%d %d %d",&tYear,&tMonth,&tDay);
    ulianT=TimeToUlian(tYear,tMonth,tDay);
    //здесь можно зациклить для вывода таблицы
    {   GelioEarth=Earth_input(ulianT,Earth);
        GeoSun.x=-GelioEarth.x;GeoSun.y=-GelioEarth.y;GeoSun.z=-GelioEarth.z;
        CoordPlanetSunEklipt=FindCoordSunEklipt(planet,ulianT);
        CoordPlanetSunEkvat=FindCoordSunEkvat(CoordPlanetSunEklipt);
        CoordEarth=FindCoordEarth(CoordPlanetSunEkvat,GeoSun);
        alpha=atan(CoordEarth.y/CoordEarth.x);
        sigma=atan(CoordEarth.z/sqrt(CoordEarth.x*CoordEarth.x+CoordEarth.y*CoordEarth.y));
        alpha=alpha/PI*180;
        sigma=sigma/PI*180;
        if(alpha<0) alpha+=360;
        AlphaClock(alpha,fd1);
        SigmaGrad(sigma,fd1);
        fclose(fd);
        fclose(fd1);
    }
    return 0;
}
