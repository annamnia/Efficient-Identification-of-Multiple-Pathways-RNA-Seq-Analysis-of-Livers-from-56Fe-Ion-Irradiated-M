#include "uni2.h"

//modularity maximization algorithm which combines power methods and Kernighan-Lin algorithm
//with fine tuning, final tuning and agglomeration step
//refer to: https://iopscience.iop.org/article/10.1088/1742-5468/2015/02/P02003/meta
comstruct uni(double **a, int n) {
    int i, j, l, flag, ng,ngtemp,gstemp, nn, np;
	int *s,*gs,*Gitemp,**S;
	double max1,max2,dQ;
	double **B, *k;
	comstruct result;

	result.community = (int *)malloc(sizeof(int) * n);
	for (i = 0; i < n; i++) result.community[i] = 0;
	result.modularity = 0;

	srand((long int) time(NULL));

	G=(int **)malloc(sizeof(int *)*n);
	gs=(int *)malloc(sizeof(int)*n);
	s=(int *)malloc(sizeof(int)*n);
	Gitemp=(int *)malloc(sizeof(int)*n);
	S=(int **)malloc(sizeof(int *)*n);

	k=(double *)malloc(sizeof(double)*n);
	B=(double **)malloc(sizeof(double *)*n);
	for (i=0; i<n; i++) {
		B[i]=(double *)malloc(sizeof(double)*n);
		G[i]=(int *)malloc(sizeof(int)*n);
		S[i]=(int *)malloc(sizeof(int)*n);
	}


	m=0;
	Q=0;
	for (i=0; i<n; i++) {
		k[i]=0;
		for (j=0; j<n; j++)
			k[i]=k[i]+a[i][j];
		
		m=m+k[i];
		//printf("%f\n", k[i]);
	}

	m=m/2;
	//printf("m is %f\n\n",m);
	if (m==0) {
		printf("a is 0 Matrix!!");
		return result;
	}

	flag=0;
	for (i=0; i<n; i++)
		for (j=0; j<n; j++) {
			//S[i][j]=0;//new...................
			B[i][j]=a[i][j]-k[i]*k[j]/(2*m);
			if (B[i][j]!=0) flag=1;
		}

	if (flag==0) {
		printf("B is 0 Matrix!");
		return result;
	}

	max1=-2;
	max2=-1;

	for (i=0; i<n; i++) {
		gs[i]=0;
		Gitemp[i]=0;
		G[0][i]=i;
	}

	for (i=1; i<n; i++)
		for (j=0; j<n; j++)
			G[i][j]=0;
	nn=np=0;
	ng=1;
	gs[0]=n;


	while (max2-max1>toler) {
		max1=max2;
		ngtemp=ng;

		for (i=0; i<ngtemp; i++){
			dQ=bisection(B,G[i],s,gs[i]);
			//printf("dQ of bisec is %f\n", dQ);
			if (dQ>toler) {
				ng++;

				gstemp=gs[i];

				for (j=0; j<gstemp; j++) {
					Gitemp[j]=G[i][j];
					G[i][j]=0;
				}


				for (j=0; j<gstemp; j++) {
					if (s[j]==1) {
						G[i][np]=Gitemp[j];
						np++;//put np to 0 later
					}
					else{
						G[ng-1][nn]=Gitemp[j];
						nn++;//put nn to 0 late
					}
				}

				gs[i]=np;
				gs[ng-1]=nn;
				np=0;
				nn=0;

				for (j=0; j<gstemp; j++){
					s[i]=0;
					Gitemp[j]=0;
				}

				Q=Q+dQ;
				/*
				 printf("\nNow there are %d groups:\n",ng);
				 for (j=0; j<ng; j++) {
				 printf("Group NO.%d: ",j+1);
				 for (l=0; l<gs[j]; l++)
				 printf("%d ",G[j][l]);
				 printf(",it has %d elements.",gs[j]);
				 printf("\n");
				 }*/

				//printf("The Modularity Q is %.6f.\n\n",Q);

			}
		}

		for (i=0; i<n; i++)
			for (j=i; j<n; j++)
				S[i][j]=S[j][i]=0;
		trans(S,ng,gs,G);
        //printf("bisec, Q is %f\n", Q);
		Q=Q+Finaltuning(B,S,n,&ng);
		invertrans(S,n,&ng,gs,G);
		Q=Q+Agglomeration(B,n,&ng,gs,G);
    //printf("agg\n" );
		//printf(".........ng is %d\n\n",ng);
		max2=Q;
	}
	//printf("The Max Modularity Q is %.10f.\n",Q);
	result.modularity = Q;
	for (i=0; i<ng; i++)
		for (j=0; j<gs[i]; j++)
			result.community[G[i][j]]=i + 1;


	for (i=0; i<n; i++) {
		free(G[i]);
		free(B[i]);
		free(S[i]);
	}
	free(gs);
	free(s);
	free(Gitemp);
	free(k);



    return result;
}

double bisection(double **B, int *g, int *s, int n){
	int i,j,nn=0,np=0,flag=0;
	double x0,x,bik,dQ;
	double *V,**c;

	if (n<=1) return 0;
	V=(double *)malloc(sizeof(double)*n);
	c=(double **)malloc(sizeof(double *)*n);
	for (i=0; i<n; i++)
		c[i]=(double *)malloc(sizeof(double)*n);

	for (i=0; i<n; i++)
		for (j=i+1; j<n; j++) {
			c[i][j]=c[j][i]=B[g[i]][g[j]];
			if (c[i][j]!=0) flag=1;
		}
	for (i=0; i<n; i++) {
		bik=0;
		for (j=0; j<n; j++)
			bik=bik+B[g[i]][g[j]];
		c[i][i]=B[g[i]][g[i]]-bik;
		if (c[i][i]!=0) flag=1;
	}

	/*for (i=0; i<n; i++) {
	 for (j=0; j<n; j++)
	 printf("%.2f",c[i][j]);
	 printf("\n");
	 }*/

	if (flag==0) {
		for (i=0; i<n; i++) free(c[i]);
		free(c);
		free(V);
		return 0;
	}

    if (n == 2) {
    	s[0] = s[1] = 1;
    } else {
		x0=eigen(c, V, n);
		if (x0<0) {
			for (i=0; i<n; i++) c[i][i]=c[i][i]-x0;
			x=eigen(c, V, n)+x0;
			for (i=0; i<n; i++) c[i][i]=c[i][i]+x0;
		}
		else x=x0;

		for (i=0; i<n; i++) {
			if (V[i]>0.0) {
				np++;
				s[i]=1;
			}
			else if (V[i]<0.0) {
				nn++;
				s[i]=-1;
			}
			else if (smprng()>0.5) {
				np++;
				s[i]=1;
			}
			else {
				nn++;
				s[i]=-1;
			}
		}
	}


	dQ=Modularity(n,s,c)/(4*m);
	//printf("dQ is %.10f\n",dQ);
	dQ=dQ+KLtuning(c,s,n);
	//printf("dQ is %.10f\n",dQ);

	for (i=0; i<n; i++)
		free(c[i]);
	free(c);
	free(V);

	return dQ;
}

double eigen(double **B, double *v, int n) {
	int i,j,N=0;
	double w1=0,w2=10,sum;
	double *u,*ad,*temp;

	u=(double *)malloc(sizeof(double)*n);
	ad=(double *)malloc(sizeof(double)*n);
	temp=(double *)malloc(sizeof(double)*n);

	for (i=0; i<n; i++)
		ad[i]=B[i][i];//to preserve the diagonal elements of **a

	for (i=0; i<(int)(n/2); i++) {//the initial u and v are random
		u[i]=1;
		v[i]=1;
	}
	for (i=(int)(n/2); i<n; i++) {
		u[i]=2;
		v[i]=2;
	}
	for (i=0; i<n; i++) temp[i]=v[i]-u[i];
	w1=abmax(u,n);

	while (fabs(w1-w2)>eps|| fabs(abmax(temp,n))>eps) {//|| fabs(abmax(temp,n))>eps
		N++;
		if (N>1000) {
			for (i=0; i<n; i++)
				B[i][i]=B[i][i]+alpha;
			N=0;
		}
		for (i=0; i<n; i++) {
			sum=0;
			for (j=0; j<n; j++)
				sum+=B[i][j]*v[j];
			u[i]=sum;
		}
		w1=w2;
		w2=abmax(u, n);
		while (w2==0) {
			for (i=0; i<n; i++)
				v[i]=i+0.5;//change
			for (i=0; i<n; i++) {
				sum=0;
				for (j=0; j<n; j++)
					sum+=B[i][j]*v[j];
				u[i]=sum;
			}
			w2=abmax(u, n);
		}

		//printf("N is %d, w1 is %.2f,w2 is %.2f\n\n",N,w1,w2);
		//printf("v is :");
		//for (i=0; i<n; i++) printf("%.2g ",v[i]);
		//printf("\nu is :");

		for (i=0; i<n; i++) {
			u[i]=u[i]/w2;
			//printf("%.2g ",u[i]);
			temp[i]=u[i]-v[i];
			v[i]=u[i];
		}

		//printf("temp is:");
		//for (i=0; i<n; i++)  printf("%.2g ",temp[i]);
		//printf("\n\n");
	}
	//printf("N is %d \n\n",N);

	for (i=0; i<n; i++)
		B[i][i]=ad[i];

	free(u);
	free(ad);
	free(temp);
	return w2;
}

double abmax(double *v,int n){
	int i;
	double temp=v[0];
	for(i=1;i<n;i++)
		if(fabs(v[i])>fabs(temp) || (fabs(v[i])==fabs(temp) && v[i]>temp))
			temp=v[i];
	return temp;
}


double KLtuning(double **B, int *s, int n){
	int i, j, k, nn, p, r;
	int *al, *c, *d;
	double dQ,qq,max1,max2;
	double *W, *dW, *dq, *q;

	max1=-2;
	max2=0;

	al=(int *)malloc(sizeof(int)*n);
	c =(int *)malloc(sizeof(int)*n);
	d =(int *)malloc(sizeof(int)*n);

	W =(double *)malloc(sizeof(double)*n);
	dW=(double *)malloc(sizeof(double)*n);
	dq=(double *)malloc(sizeof(double)*n);
	q =(double *)malloc(sizeof(double)*n);

	while (max2-max1>toler) {
		//for (i=0; i<n; i++)
		//	printf("%d ",s[i]);
		//printf("\nmax1 is %.3f, max2 is %.3f\n",max1,max2);
		max1=max2;
		q[0]=max2;

		for (i=0; i<n; i++) {
			al[i]=-1;
			c[i]=0;//c[i]=0 means node i not moved yet
			dW[i]=0;
			W[i]=0;
		}

		for (i=0; i<n; i++) {
			for (j=0; j<n; j++)
				W[i]=W[i]+(double)s[j]*B[i][j];
			//printf("W is %.10f\n",W[i]);
			dq[i]=(B[i][i]-(double)s[i]*W[i])/m;//calculate dQ
			//printf("for node %d,dq is %.10f\n",i,dq[i]);
		}
		//printf("\n\n");

		/* to pick up the largest dQ*/
		dQ=-10;
		for (i=0; i<n; i++) {
			if (dq[i]>dQ) {
				p=1;
				d[0]=i;
				dQ=dq[i];
			}
			else if (dq[i]==dQ) {
			    d[p]=i;
				p++;
			}
		}

		if (p==1) {
			al[0]=d[0];
			c[al[0]]=1;//c[i]=1 means the i-th node has been choosen
		}
		else {
			p=(int)p*smprng();
			al[0]=d[p];
			c[al[0]]=1;
		}
		//s[al[0]]=-s[al[0]];
		//printf("NO.0 dQ is %.10f,node %d is moved\n",dQ,al[0]);
		//for (i=0; i<n; i++)
		//	printf("%d ",s[i]);
		//printf("\n\n");

		/*moving nodes*/
		for (i=1; i<n; i++) {
			q[i]=q[i-1]+dQ;
			dQ=-10;
			for (j=0; j<n; j++)
				if (c[j]==0) {
					dW[j]=(double)s[al[i-1]]*B[j][al[i-1]]*2;
					//printf("dW is %.10f\n",dW[j]);
					W[j]=W[j]-dW[j];
					dq[j]=(B[j][j]-(double)s[j]*W[j])/m;
					//printf("for node %d,dq is %.10f\n",j,dq[j]);

					if (dq[j]>dQ) {
						p=1;
						d[0]=j;
						dQ=dq[j];
					}
					else if (dq[j]==dQ) {
						d[p]=j;
						p++;
					}
				}
			//printf("\n\n");
			/*
			for (j=0; j<n; j++)
				if (c[j]==0) {
					if (dq[j]>dQ) {
						p=1;
						d[0]=j;
						dQ=dq[j];
					}
					else if (dq[j]==dQ) {
						d[p]=j;
						p++;
					}
				}*/
			//printf("dQ is %.10f\n",dQ);
			if (p==1) {
				al[i]=d[0];
			}
			else {
				p=(int)p*smprng();
				al[i]=d[p];
			}
			c[al[i]]=1;
			//s[al[i]]=-s[al[i]];
			//printf("NO.%d dQ is %.10f,node %d is moved\n",i,dQ,al[i]);
			//for (j=0; j<n; j++)
			//	printf("%d ",s[j]);
			//printf("\n\n");
		}

		/*to find the largest configuration*/
		qq=-1;
		for (i=0; i<n; i++) {
			//printf("Q is %.3f\n",q[i]);
			if (q[i]>qq) {
				p=1;
				d[0]=i;
				qq=q[i];
			}
			else if (q[i]==qq) {
				d[p]=i;
				p++;
			}
		}
		if (p==1)
			r=d[0];
		else {
			p=(int)p*smprng();
			r=d[p];
		}//the #r-th state has the largest configuration
		//printf("r is %d\n\n\n",r);
		for (i=0; i<r; i++)
			s[al[i]]=-s[al[i]];
		dQ=q[r]-max2;
		max2=q[r];
		//printf("one round of dQ is %.10f\n",dQ);
	}

	free(al);
	free(c);
	free(d);
	free(W);
	free(dW);
	free(dq);
	free(q);

	return max1;
}

double Modularity(int n,int *s,double **b){
	int i,j;
	double q=0;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			q=q+(double)s[i]*b[i][j]*s[j];

	return(q);
}



double Addmatrix(double **b,int n){
	int i,j;
	double sum=0;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			sum=sum+b[i][j];

	return(sum);
}

void trans(int **S, int ng, int *gs, int **G){
	int i,j;

	for (i=0; i<ng; i++)
		for (j=0; j<gs[i]; j++)
			S[G[i][j]][i]=1;

	return;
}

void invertrans(int **S, int n, int *ng, int *gs, int **G){
	int i,j,k,temp,flag=0;

	temp=*ng;

	for (i=0; i<n; i++) {
		gs[i]=0;
		for (j=0; j<n; j++)
			G[i][j]=0;
	}//initialize gs and G

	k=0;
	for (i=0; i<*ng; i++) {
		for (j=0; j<n; j++)
			if (S[j][i]==1) {
				G[k][gs[k]]=j;
				gs[k]++;
				flag=1;
			}
		if (flag==1) k++;
		flag=0;
	}

	//printf("k is %d \n",k);
	*ng=k;

	return;

}

double Finaltuning(double **B, int **S, int n, int *ng){
	int i,j,k, p, r;
	int *con, *al, *be, *c, *dn, *dg;
	double dQ, qq, max1, max2;
	double **W, **dW, **dq, *q;

	max1=-2;
	max2=0;

	con=(int *)malloc(sizeof(int)*n);//one dimensional configuration
	al=(int *)malloc(sizeof(int)*n);//record the node's to be moved each time
	be=(int *)malloc(sizeof(int)*n);//record the group each node's being moved to
	c =(int *)malloc(sizeof(int)*n);
	dn=(int *)malloc(sizeof(int)*n*n);
	dg=(int *)malloc(sizeof(int)*n*n);

	W=(double **)malloc(sizeof(double *)*n);
	dW=(double **)malloc(sizeof(double *)*n);
	dq=(double **)malloc(sizeof(double *)*n);
	q =(double *)malloc(sizeof(double)*n);

	for (i=0; i<n; i++) {
		W[i]=(double *)malloc(sizeof(double)*n);
		dW[i]=(double *)malloc(sizeof(double)*n);
		dq[i]=(double *)malloc(sizeof(double)*n);
	}

	for (i=0; i<n; i++) {
		j=0;
		while (S[i][j]==0) {
			j++;
		}
		con[i]=j;
	}

	while (max2-max1>toler) {
		//printf("\nmax1 is %.3f, max2 is %.3f\n",max1,max2);
		max1=max2;
		q[0]=max2;

		for (i=0; i<n; i++) {
			al[i]=-1;
			c[i]=0;
			for (j=0; j<n; j++) {
				dW[i][j]=0;
				W[i][j]=0;
			}
		}

		for (i=0; i<n; i++)
			for (j=0; j<*ng+1; j++)
				for (k=0; k<n; k++)
					W[i][j]=W[i][j]+(double)B[i][k]*S[k][j];

		for (i=0; i<n; i++)
			for (j=0; j<*ng+1; j++) {
				if (j==con[i]) dq[i][j]=0;
				else {
					dq[i][j]=(B[i][i]-W[i][con[i]]+W[i][j])/m;
				//    printf("for node %d to group %d ,dq is %.10f\n",i,j,dq[i][j]);
				}
			}

		dQ=-10;
		for (i=0; i<n; i++)
			for (j=0; j<*ng+1; j++)
				if (j!=con[i]) {
					if (dq[i][j]>dQ) {
						p=1;
						dn[0]=i;
						dg[0]=j;
						dQ=dq[i][j];
					}
					else if (dq[i][j]==dQ) {
						dn[p]=i;
						dg[p]=j;
						p++;
					}
				}

		if (p==1) {
			al[0]=dn[0];
			c[al[0]]=1;
			be[0]=dg[0];
		}
		else {
			p=(int)p*smprng();
			al[0]=dn[p];
			c[al[0]]=1;
			be[0]=dg[p];
		}

		if (be[0]==*ng) *ng=*ng+1;

		//printf("NO.0 dQ is %.10f,node %d is moved to group %d\n\n\n",dQ,al[0],be[0]);

		/*moving nodes*/
		for (i=1; i<n; i++) {
			q[i]=q[i-1]+dQ;
			dQ=-10;

			for (j=0; j<n; j++)
				if (c[j]==0)
					for (k=0; k<*ng+1; k++) {
						if (k==con[al[i-1]]) dW[j][k]=dW[j][k]+B[j][al[i-1]];
						else if (k==be[i-1]) dW[j][k]=dW[j][k]-B[j][al[i-1]];
						//printf("dW is %.10f\n",dW[j][k]);
						//W[j][k]=W[j][k]-dW[j][k];
					}

			for (j=0; j<n; j++)
				if (c[j]==0)
					for (k=0; k<*ng+1; k++) {
						if (k==con[j]) dq[j][k]=0;
						else {
							dq[j][k]=(B[j][j]-W[j][con[j]]+W[j][k]+dW[j][con[j]]-dW[j][k])/m;
						//	printf("for node %d to group %d,dq is %.10f\n",j,k,dq[j][k]);
						}
					}

			for (j=0; j<n; j++)
				if (c[j]==0)
					for (k=0; k<*ng+1; k++)
						if (k!=con[j]) {
							if (dq[j][k]>dQ) {
								p=1;
								dn[0]=j;
								dg[0]=k;
								dQ=dq[j][k];
							}
							else if (dq[j][k]==dQ) {
								dn[p]=j;
								dg[p]=k;
								p++;
							}
						}

			if (p==1) {
				al[i]=dn[0];
				be[i]=dg[0];
			}
			else {
				p=(int)p*smprng();
				al[i]=dn[p];
				be[i]=dg[p];
			}
			c[al[i]]=1;

			if (be[i]==*ng) *ng=*ng+1;
			//printf("NO.%d dQ is %.10f,node %d is moved to group %d\n\n\n",i,dQ,al[i],be[i]);
		//	printf("group number is %d\n\n",*ng);
		}



		qq=-1;
		for (i=0; i<n; i++) {
		//	printf("Q is %.3f\n",q[i]);
			if (q[i]>qq) {
				p=1;
				dn[0]=i;
				qq=q[i];
			}
			else if (q[i]==qq) {
				dn[p]=i;
				p++;
			}
		}

		if (p==1) r=dn[0];
		else {
			p=(int)p*smprng();
			r=dn[p];
		}
		//printf("r is %d\n\n\n",r);
		for (i=0; i<r; i++){
			S[al[i]][con[al[i]]]=0;
			S[al[i]][be[i]]=1;
            con[al[i]]=be[i];
		}

		dQ=q[r]-max2;
		max2=q[r];
		//printf("one round of dQ is %.10f\n",dQ);
	}

	free(al);
	free(be);
	free(c);
	free(dn);
	free(dg);
	free(con);
	free(q);

	for (i=0; i<n; i++) {
		free(W[i]);
		free(dW[i]);
		free(dq[i]);
	}
	free(W);
	free(dW);
	free(dq);

	return max1;
}


double Agglomeration(double **B, int n, int *ng, int *gs, int **G) {
	int i,j,k,R,commnum,g1,g2;
	double dQ,temp,max;
	int *con, **config;
	double *Qt;

	commnum=*ng;
	Qt=(double *)malloc(sizeof(double)*n*n/2);
	con=(int *)malloc(sizeof(int)*n);
	config=(int **)malloc(sizeof(int *)*commnum);
	for (i=0; i<commnum; i++) config[i]=(int *)malloc(sizeof(int)*n);

	Qt[0]=0;

	for (i=0; i<commnum-1; i++) {
		/*record */
		for (j=0; j<*ng; j++)
			for (k=0; k<gs[j]; k++)
				config[i][G[j][k]]=j;
		for (j=0; j<n; j++) con[j]=config[i][j];

		temp=-1;
		for (j=0; j<*ng-1;j++ )
			for (k=j+1; k<*ng; k++) {
				dQ=dQA(j,k,gs,G,B);
				if (dQ>temp) {
					temp=dQ;
					g1=j;
					g2=k;
				}
			}

		Qt[i+1]=Qt[i]+temp;
		for (j=0; j<n; j++) {
			if (con[j]==g2)
				con[j]=g1;
			else if (con[j]>g2)
				con[j]=con[j]-1;
		}
		con_to_G(con,n,ng,gs,G);
	}

	for (i=0; i<*ng; i++)
		for (j=0; j<gs[i]; j++)
			config[commnum-1][G[i][j]]=i;

	temp=-1;
	for (i=0; i<commnum; i++)
		if (Qt[i]>=temp) {
			temp=Qt[i];
			R=i;
		}

	con_to_G(config[R],n,ng,gs,G);
	max=Qt[R];

	for (i=0; i<commnum; i++)
		free(config[i]);
	free(config);
	free(Qt);

	return max;
}

/*dQ in the agglomeration step*/
double dQA(int n1, int n2, int *g1, int **G, double **b) {
	int i,j;
	double dQ=0;

	for (i=0; i<g1[n1]; i++)
		for (j=0; j<g1[n2]; j++)
			dQ=dQ+b[G[n1][i]][G[n2][j]];
	dQ=dQ/m;

	return dQ;
}


void con_to_G(int *con, int n, int *ng, int *gs, int **G) {
	int i,j;

	*ng=0;
	for (i=0; i<n; i++) {
		gs[i]=0;
		for (j=0; j<n; j++)
			G[i][j]=0;
	}//initialize gs and G

	for (i=0; i<n; i++) {
		G[con[i]][gs[con[i]]]=i;
		gs[con[i]]++;
		if (con[i]>*ng-1)
			*ng=con[i]+1;
	}

	return;
}



double smprng(void){

	return ((double) rand())/((double) RAND_MAX);

}
