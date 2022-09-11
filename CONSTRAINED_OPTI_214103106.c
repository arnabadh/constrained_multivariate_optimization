#include<stdio.h>
#include<conio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define epsilon pow(10, -3)

double Bounding_Phase(double a[100], double S[100], int n);
double objective_function(double a[100],double S[100],double x, int n);
double func(double x[100],int n);
double g_x(double x[100],int n);
double R=0.1;
double Golden_Section_Search(double z[100]/*the x point*/,double S[100]/*search direction*/,double a, double b,int n /*number of variable*/);
double penalty(double g);


double objective_function(double a[100]/*the starting point*/,double S[100],double x, int n)
{ 
int i;
double c[n+1];
for(i=1;i<=n;i++)
{
    c[i]=a[i]+x*S[i];
}

	return func(c,n);
}

double Bounding_Phase(double a[100] , double S[100], int n)/*the starting point*/
{
     
	int i=1 , feval,j=0;/*Number of iterations*/
	double delta=0.1;
	double x0=0.2, f_1, f0, f1,x[50000],f[50000];
	 
	/*Step 1*/
	feval = 0; /*function evaluation*/
	f_1 = objective_function(a,S,x0-delta,n); /*Calculate objective_function*/
	f0 = objective_function(a,S,x0,n);
	f1 = objective_function(a,S,x0+delta,n);
	
	if(f_1>f0 && f0>f1) 
	delta=delta;
	else if(f_1<f0 && f0<f1) 
	delta=-delta;
	x[0]=x0;
	f[0]=f0;

	 
	feval = feval + 3;
	

	do{  j=j+1;
		
		x[j]=x[j-1]+pow(2,j-1)*delta;
		f[j]=objective_function(a,S,x[j],n);
		
		if(f[j]>f[j-1])
		break;
		feval++;
		i = i + 1;
		  
	
	}while( (f[j]<f[j-1] || f[j]==f[j-1]) && i <20);/*Step 3*/

	
	return Golden_Section_Search(a,S,x[j-2],x[j],n);//calling golden search section with average of value got in bounding phase
	
}


double Golden_Section_Search(double z[100]/*the x point*/,double S[100]/*search direction*/,double a, double b,int n /*number of variable*/)
{
  double aw=0,bw=1,lw=1,e=0.00001,fw1,fw2,w1=0,w2=0,x1=0,x2=2,xavg,fxavg,c;
  int i=1,feval=1;
 //printf("\na = %lf  -   b = %lf\n",a,b);
 if(a>b){
 	c=a;
 	a=b;
 	b=c;
 }
 //printf("3");
  do{
      
      w1=aw+0.618*lw;//normalising the function with respect to w
      w2=bw-0.618*lw;
      
      x1=w1*(b-a)+a;
      x2=w2*(b-a)+a;
   
      fw1=objective_function(z,S,x1,n);
      fw2=objective_function(z,S,x2,n);
		
      if(fw2<fw1)//region elimination algo
          {
              bw=w1;
              lw= fabs(bw-aw);
          }
      else if(fw1<fw2)
          {
              aw=w2;
              lw= fabs(bw-aw);
          }
          else
          { bw=w1;
            aw=w2;
            lw= fabs(bw-aw);  
          }
      
      feval++;
      i++;
      x1=aw*(b-a)+a;
      x2=bw*(b-a)+a;
  }while(lw>e );//termination condition of golden search section
 //printf("4");
 	
    return (x1+x2)/2.0;//returns the minimum value of alpha
  
}


double *powell_conjugate_direction(double x[200], int n)
{
	
	double  x_updated[200],alpha,d[200],dabs=0,t=0,* x_pointer;
    double S[200][200];
    int i,j,Ns=0;
    
    x_pointer = malloc(sizeof(double)*(n+2));
    
    for(i=1;i<=n;i++)//intialising the n orthogonal search direction
    {
    for(j=1;j<=n;j++)
    {
    	
    if(i==j)
    S[j][i]=1;
    
    else
    S[j][i]=0;
    
    }
    }
    
    
    
    double S_old[n+1], store[n+1][3];//S_old to store the search direction in 1D array //store to store the value of recent x point 
    
    for(i=1;i<=n;i++)
    x_updated[i]=x[i];
    
   // fprintf(op1,"\n%d\t%lf",Ns,func(x_updated, n));//printing the number of iteration vs function value in output file
    
	do{
		//printf("\n1");
		Ns++;//counter for iteration
		
	   dabs=0;
	   
       for(i=1;i<=n;i++)
       {  
        S_old[i]=S[i][1];//storing the search direction along FIRST SEARCH DIRECTION in S_old
       }
   
    	
    
        alpha= Bounding_Phase(x_updated,S_old,n);//calculating the value of alpha
       // printf("\n%lf",alpha);//printing the value of alpha in console
      //  fprintf(op2,"\n%lf",alpha);//printing the value of alpha in output file
        
        //printf("\n2");
        for(i=1;i<=n;i++){
        x_updated[i]= x_updated[i]+ alpha*S_old[i];//evaluating value of next x corresponding to given alpha
        store[i][1]=x_updated[i];
     //   printf("\t%lf",x_updated[i]);//printing the value of x in console
     //   fprintf(op2,"\t%lf",x_updated[i]);//printing the value of x in output file
        }
        
       
        
        
    
       // printf("\n");
     //   fprintf(op2,"\n");
        
    
    
        for(j=2;j<=n;j++)//calculating ALPHA AND CORRESPONDING X IN CORRESPONDING SEARCH DIRECTIONS
        {
    	
         for(i=1;i<=n;i++)
		 {
		 S_old[i]=S[i][j];
         store[i][2]=x_updated[i];
            }
         alpha= Bounding_Phase(x_updated,S_old,n);
    
    
      //  printf("\n%lf",alpha);
     //    fprintf(op2,"\n%lf",alpha);
         
         for(i=1;i<=n;i++)
		 {
        x_updated[i]= store[i][2]+ alpha*S_old[i];
     //  printf("\t%lf",x_updated[i]);
     //   fprintf(op2,"\t%lf",x_updated[i]);
        }
        
        
        
        
        }
        
        for(i=1;i<=n;i++)
        store[i][2]=x_updated[i];
        
        for(i=1;i<=n;i++)
            S_old[i]=S[i][1];//storing the search direction along FIRST SEARCH DIRECTION in S_old
        alpha= Bounding_Phase(x_updated,S_old,n);//calculating the value of alpha
     //   printf("\n%lf",alpha);//printing the value of alpha in console
     //   fprintf(op2,"\n%lf",alpha);//printing the value of alpha in output file
        
        
        
        for(i=1;i<=n;i++){
        x_updated[i]= store[i][2]+ alpha*S_old[i];//evaluating value of next x corresponding to given alpha
     //   printf("\t%lf",x_updated[i]);//printing the value of x in console
     //   fprintf(op2,"\t%lf",x_updated[i]);//printing the value of x in output file
        }
        
     //   fprintf(op1,"\n%d\t%lf",Ns,func(x_updated, n));//printing the number of iteration vs function value in output file
 
    
	    for(i=1;i<=n;i++)
	    d[i]= (x_updated[i]-store[i][1]);
	   
	    for(i=1;i<=n;i++)
	    t = t + pow( (x_updated[i]-store[i][1]), 2);
	   
	    dabs=sqrt(t);//FINDING NORM OF difference between Xn+1 and X1
	   
	    for(j=n;j>=2;j--)   //updating corresponding search directions
	    {
	    for(i=1;i<=n;i++)
	    S[i][j]=S[i][j-1];
	    }
	   
	    for(i=1;i<=n;i++)
	    S[i][1]=d[i]/dabs;//updating S1 direction
    
    	
   
   
}while(dabs>epsilon && Ns<300);//termination condition of powell conjugate
  
  
   for(i=1;i<=n;i++){
   
   x_pointer[i]=x_updated[i];
   //printf("\n%lf",x_updated[i]);
}
   return x_pointer;  
}

double func(double x[100],int n){ /*function defination of multi variable function*/
    int i,j;
    double func_value=0;
   // func_value=pow(x[1]-10,3)+pow(x[2]-20,3);
    
  //  func_value=-1*(pow(sin(6.283185307179586*x[1]),3)*sin(6.283185307179586*x[2]))/(pow(x[1],3)*(x[1]+x[2]));
    
    func_value=x[1]+x[2]+x[3];
    
    return (func_value+ R* g_x(x,n) );
    
    
}

double g_x(double x[100],int n){
	double g1,g2,g3,g4,g5,g6;
	
	
	
   //g1=pow(x[1]-5,2)+pow(x[2]-5,2)-100;
   // g2=-1*pow(x[1]-6,2)-1*pow(x[2]-5,2)+82.81;
 	
 	
 //	g1=-x[1]*x[1]+x[2]-1;
//	g2=-1+x[1]-(x[2]-4)*(x[2]-4);
	
	
	g1=-1+0.0025*(x[4]+x[6]);
	
    g2=-1+0.0025*(-1*x[4]+x[5]+x[7]);

    g3=-1+0.01*(-1*x[6]+x[8]);

    g4=100*x[1]-x[1]*x[6]+833.33252*x[4]-83333.333;

    g5=x[2]*x[4]-x[2]*x[7]-1250*x[4]+1250*x[5];

    g6=x[3]*x[5]-x[3]*x[8]-2500*x[5]+1250000;
    
    g1=-1*g1;
    g2=-1*g2;
    g3=-1*g3;
    g4=-1*g4;
    g5=-1*g5;
    g6=-1*g6;
   	return penalty(g1)+penalty(g2)+penalty(g3)+penalty(g4)+penalty(g5)+penalty(g6);
   	
   	
   	return penalty(g1)+penalty(g2);
}

double penalty(double g)
{	
	if (g<0)
	g=pow(g,2);
	
	else g=0;
	
	return g;
}



int main()
{
    double x[200],* x_pointer,x_new[200],f_r_old,f_r_new,R1;
     FILE *op;
    op = fopen("convergence_plot.txt","w");
    fprintf(op,"#It\t f_value ");
    int i,n,j,k=0;
    
    printf("Enter the number of variables ");
    scanf("%d",&n);//entering the number of variables
    
     printf("\nEnter the starting point ");
    for(i=1;i<=n;i++)
    scanf("%lf",&x[i]);//entering the starting point
    
     

    x_pointer=malloc(sizeof(double)*(n+2));
   	x_pointer= powell_conjugate_direction(x,n);
    
    
    
    printf("\n");
	for(i=1;i<=n;i++){
	   x_new[i]=x_pointer[i];
	   printf("%lf\t",x_new[i]);
	}
   
   printf("\n");
   f_r_old=func(x_new,n);
   R=R*10;
   f_r_new=func(x_new,n);
  
   
  
   while(fabs(f_r_old-f_r_new)>=0.001 && k<=50){
   	x_pointer= powell_conjugate_direction(x_new,n);
   	k=k+1;
    for(i=1;i<=n;i++){
   		x_new[i]= x_pointer[i];
      	printf("%lf\t",x_new[i]);
   	}
   printf("\n");
   f_r_new=func(x_new,n);
   R=R*10;
   f_r_old=f_r_new;
   f_r_new=func(x_new,n);
   
   R1=R;
   R=0;
   fprintf(op,"\n%d\t%lf",k,func(x_new,n));
   printf("%d ",k);
   R=R1;
   
   }
   
   
   
    for(i=1;i<=n;i++)
    {
    	printf("%lf    ",x_new[i]);
	}
	R=0;
	printf("\n%lf    ",func(x_new,n));
	
    fclose(op);
 // fclose(op2);
   return 0;
}
