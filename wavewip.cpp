# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <conio.h>

int main()
{   float a,x1,x2,T,h,tau,c;
	int i,j;
	short int n=200; 
	float u[n],uplust[n],uminust[n];
	
	printf("PDE d^2U/dx^2-(1/a^2)*d^2U/dt^2=f(x,t) approximates moving of the wave in 1D space.\n");
	printf("Programm solves hyperbolic PDE using finite differences method.\n");
	printf("\nNote: Functions m1(x)=U(x,0); m2(x)=U't(x,0); m3(t)=U(x1,t); m4(t)=U(x2,t); f(x,t) are defined in a source code of the programm.\n");   //m1(x)=x^2-4*x-9 ; m2(x)=2*x; m3(t)=t; m4(t)=t+1; f(x,t)=0	
	printf("Note: Space coordinate x1 is defined in a source code.\n");
	printf("Note: if Courant condition wont be appropriate, you'll have to reenter the time of the calculaion and the phase speed of the wave.\n");
	x1=1; //m4(t)-m3(t),t=0
	
	FILE *f1;
    if ((f1=fopen("D:\\MAINSTORAGE\\GUAP\\INF\\3S\\wavedata.txt", "w+"))==NULL)
    { 
	printf("Opening error\n");
    return 1;  	
    }
    else {
    printf("\nFile wavedata.txt for the solutions of PDE opened for writing.\n");
    }
	
	do {
	printf("\nEnter square of the phase speed of the wave.\n");
	scanf("%f",&a);
	printf("Enter the time of the calculation.\n");
	scanf("%f",&T);	
	printf("Enter the x2 space coordinate.\n");
	scanf("%f",&x2);	
	h=(x2-x1)/(n-1);
	tau=T/(n-1);
	c=(tau*sqrt(a))/h;
	printf("\nCourant number: %f\n",c);
	} while ((fabs(c)>=1)||(T<=0)||(x2<=1)); //устойчивость схемы, проверка условия Куранта, защита от ввода отрицательного времени
	
	for(j=0;j<n; j++)     
		{ uminust[j]=x1+j*h; // U(x,0), 0-й слой, из н.у
	      fprintf(f1,"%f ",uminust[j]);	
		}		
	u[0]=x1+tau; 
	u[n-1]=uminust[n-1]+tau; 
	fprintf(f1,"\n%f ",u[0]);
    for(j=1;j<n-1; j++) 
    {
    	u[j]=pow((x1+j*h),2)-4*(x1+j*h)-9+tau*2*(x1+j*h)+pow(tau,2)/2*a;  //d^2m1(x)/dx^2=1, u(x,tau), 1-й слой, из н.у. из производной по Тэйлору с погрешностью Облако(тау^3)
    	fprintf(f1,"%f ",u[j]);
	}
    fprintf(f1,"%f",u[n-1]);
    
	do
    {
	uplust[0]=u[0]+tau; 
	fprintf(f1,"\n%f ",uplust[0]);
	uplust[n-1]=u[n-1]+tau;
    for(j=1;j<n-1; j++)
    { uplust[j]=2*(1-c*c)*u[j]-uminust[j]+c*c*(u[j+1]+u[j-1]);
	 fprintf(f1,"%f ",uplust[j]);
	}	
	fprintf(f1,"%f",uplust[n-1]);
	
	for(j=0;j<n; j++) //перезапись вычисленного в текущий
    {
    uminust[j]=u[j];
    u[j]=uplust[j];
   }
   	} while (uplust[0]<T);	
	printf("\nAll data succesfully written to wavedata.txt.\n");
	getch();
	return 0;
}
