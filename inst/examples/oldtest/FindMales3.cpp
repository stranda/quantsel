#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

#define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))


using namespace Rcpp;
using namespace sugar;


// [[Rcpp::export]]
IntegerVector oneMultinomCall(NumericVector probs) {
    int k = probs.size();
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    return(ans);
}


// [[Rcpp::export]]
LogicalVector matchl(NumericVector X, NumericVector i) {
  NumericVector m(X.size());
  LogicalVector l(X.size());
  m=match(X,i);
  for (int j=0;j<m.size();j++) {if (m(j)>0) l(j)=TRUE; else l(j)=FALSE;}
  return l;
}


// rcpp in; rcpp out
// [[Rcpp::export]]
NumericMatrix submat(NumericMatrix X, LogicalVector condition) { 
    int n=X.nrow(), k=X.ncol();
    NumericMatrix out(sum(condition),k);
    for (int i = 0, j = 0; i < n; i++) {
        if(condition[i]) {
            out(j,_) = X(i,_);
            j = j+1;
        }
    }
    return(out);
}


// [[Rcpp::export]]
LogicalVector InCircle_sugar(NumericVector circle, NumericMatrix points)
{ //slower than the 'for loop approach' below
  double cx=circle[0],cy=circle[1],cr=circle[2];
  double cr2 = pow(cr,2);

  //uses the vectorized nature of Rcpp sugar additions
  
  return   ( (points(_,0)-cx)*(points(_,0)-cx) - (points(_,1)-cy)*(points(_,1)-cy) <= cr2 );
}



// [[Rcpp::export]]
LogicalVector InCircle1(NumericVector circle, NumericVector x, NumericVector y)
{
  double cx=circle[0],cy=circle[1],cr=circle[2];
  double cr2 = pow(cr,2);
  return (((x-cx)*(x-cx) + (y-cy)*(y-cy))<=cr2);
}

// [[Rcpp::export]]
LogicalVector InCircle(NumericVector circle, NumericMatrix points)
{
  LogicalVector ret(points.nrow());
  double cx=circle[0],cy=circle[1],cr=circle[2];
  int i ;
  double cr2 = pow(cr,2);

  for (i=0;i<points.nrow();i++)
    {
      ret[i]=(((points(i,0)-cx)*(points(i,0)-cx) + (points(i,1)-cy) * (points(i,1)-cy))<=cr2);
    }
  
  return ret;
}



// [[Rcpp::export]]
NumericMatrix FindMales3(NumericVector reprows,
			 NumericMatrix individuals,
			 NumericVector mrdist,
			 NumericVector tries,
			 NumericVector pct,
			 NumericVector deltapct)

{
  int cnt = 0;
  int dadcnt = 0;
  int complete = 0;
  int i, j, r, rint;
  double d, tmpdst;
  NumericMatrix ret(reprows.size(),3);
  ret(_,0) = reprows;
  
  LogicalVector icircO(individuals.nrow());
  LogicalVector icircI(individuals.nrow());
  LogicalVector completed(reprows.size());
  
  for (i=0;i<completed.size();i++) completed[i]=FALSE;
  
  NumericVector circO(3);
  NumericVector circI(3);

  NumericVector dads(10000);

  
  
  

  NumericMatrix points (individuals.nrow(),2) ;

  points(_,0) = individuals(_,3);
  points(_,1) = individuals(_,4);
  

  while ((cnt<tries(0)) & (complete==0))
    {
      for (i=0; i<reprows.size(); i++)
	{
	  if (!completed[i])
	    {

	      r = reprows[i]-1;
	  
	      circO[0]=individuals(r,3);
	      circO[1]=individuals(r,4);
	      circO[2]=mrdist(i) * (1 + pct(0));
	      circI[0]=circO[0];
	      circI[1]=circO[1];
	      circO[2]=mrdist(i) * (1 - pct(0));

	      icircO = InCircle(circO,points);
	      icircI = InCircle(circI,points);

	      dadcnt = 0;
	      for (j=0; j<points.nrow();j++)
		{
		  if ((icircO[j])&(!icircI[j]))
		    {
		      dads[dadcnt]=j;
		      dadcnt++;
		    }
		}
	      if (dadcnt>0)
		{
		  rint = 0;
		  d = 10000000;
		  for (j=0;j<dadcnt;j++)
		    {
		      tmpdst = sqrt(pow(points(r,0)-points(dads(j)-1,0),2)+
					   pow(points(r,1)-points(dads(j)-1,1),2));
		      if (abs(mrdist(i)-tmpdst) < d )
			{
			  d = tmpdst;
			  rint=j;
			}
		    
		    }
		  ret(i,1) = dads(rint)+1;
		  ret(i,2) = sqrt(pow(points(r,0)-points(ret(i,1)-1,0),2)+
				  pow(points(r,1)-points(ret(i,1)-1,1),2));
		  //		  Rprintf("points for mom: %g, %g and points %g, %g for dad\n",points(r,0),points(r,1),points(dads(rint),0),points(dads(rint),1));
		  completed[i] = TRUE;
		  
		}
	    }
	  // Rprintf("made it through one mothers %i search for dad. found %i dads\n",r,dadcnt);  
	}

      //      Rprintf("here\n");
      
      if (sum(!completed)==0) complete = 1;   
      
      pct(0) = pct(0) + deltapct(0);
      cnt++ ;
    }
  
    return ret;
}


// [[Rcpp::export]]
NumericVector FindMales4(IntegerVector reprows,
			 NumericMatrix individuals,
			 NumericVector radii,
			 NumericMatrix pollenkern)
  
{
  //  Timer timer;
  int i, j, k, r;
  double x,y,rd;
  
  NumericMatrix ret(reprows.size(),2);
  ret(_,0) = reprows;
  
  LogicalVector m(individuals.nrow());
  
  NumericVector circO(3);
  
  NumericVector dads(10000);
  NumericVector dsts(individuals.nrow());
  NumericVector dens(individuals.nrow());
  
    NumericMatrix points (individuals.nrow(),4) ;
  NumericMatrix points2 (individuals.nrow(),4) ;

  points2(_,0) = individuals(_,3); //x coord
  points2(_,1) = individuals(_,4); //y coord
  points2(_,2) = seq(0,individuals.nrow()-1); //keep original row numbers to use as father ids in output
  points2(_,3) = individuals(_,0); //democlass

  
  for (i=0; i<reprows.size(); i++)
    {//figure out which pops are overlaid by the circle. only keep those individuals
      r = reprows(i)-1;
      x=individuals(r,3);
      y=individuals(r,4);
      rd=radii(i);
      circO(0) = x;
      circO(1) = y;
      circO(2) = rd;
      
      m = InCircle(circO,points2);
      points=submat(points2,m) ;  
      
      if (points.nrow()>0)
	{
	  dsts = sqrt(pow(x-points(_,0),2) + pow(y-points(_,1),2));
	  if (pollenkern(1,0)==3) //mixed dist or bust
	    {
	      dens = dweibull(dsts,pollenkern(1,1),
			      pollenkern(1,2))*(1-pollenkern(1,5)) +
		dnorm(dsts,pollenkern(1,3),
		      pollenkern(1,3))*pollenkern(1,5);
	    } else Rcpp::stop("pollen kernel other than 3 specified");
	  
	  dens = dens/sum(dens); //normalize
	  
	  ret(i,1) = points(which_max(oneMultinomCall(dens)),2) + 1; //added one to return to R indexing
	}
      else
	{
	  ret(i,1)= -1 ;
	}
    }

  return ret;
}

// [[Rcpp::export]]
NumericVector FindMalesSTL(IntegerVector reprows,
			   NumericMatrix individuals,
			   NumericVector populations, //population for each individual 1-indexed
			   NumericVector radii,
			   IntegerVector stgs,
			   NumericMatrix locs,     //leftx,boty,rightx,topy
			   NumericMatrix pollenkern)

{//this works as is but is 3 or 4 times slower than the FindMales4

  //  Timer timer;
  int i, j, k, r, dad, inds;
  double x,y,rd, dst,ds;
  dst=0;
  NumericMatrix ret(reprows.size(),2);
  ret(_,0) = reprows;

  LogicalVector m(individuals.nrow());
  
  NumericVector circO(3);
  
  NumericVector dsts(individuals.nrow());
  std::vector< double > dens;
  std::vector< int > oid;

  dens.reserve(individuals.nrow());
  oid.reserve(dens.size());

  
  std::vector< std::vector< double> > points(individuals.nrow(), std::vector<double>(4));//should be points[r][c]
  NumericMatrix points2(individuals.nrow(),4);

  for (i=0;i<individuals.nrow();i++)
      {
	points[i][0] = individuals(i,3);
	points[i][1] = individuals(i,4);
	points[i][2] = i;
	points[i][3] = individuals(i,0);
      }
  inds=individuals.nrow();
 
  for (i=0; i<reprows.size(); i++)
    {//figure out which pops are overlaid by the circle. only keep those individuals

      //timer.step("start");
      
      r = reprows(i);
      x=individuals(r,3);
      y=individuals(r,4);
      rd=radii(i);
      circO(0) = x;
      circO(1) = y;
      circO(2) = rd;
      


      dens.resize(0);
      oid.resize(0);

      m=InCircle1(circO,individuals(_,3),individuals(_,4));

      for (j=0;j<inds;j++)
	{
	  //	  Rprintf("j %i \n",j);
	  if (m[j])
	    {
	      dst = sqrt(pow(x-points[j][0],2) + pow(y-points[j][1],2));
	      //	      Rprintf("j %i ->",j);
	      //	      Rprintf("dst %g \n",dst);
	      if (pollenkern(1,0)==3) //mixed dist or bust
		{
		  dens.push_back (R::dweibull(dst,pollenkern(1,2), pollenkern(1,1),0) * (1-pollenkern(1,5)) +
				  R::dnorm(dst,pollenkern(1,3),pollenkern(1,3),0)*pollenkern(1,5));
		  oid.push_back(j+1);
	      } else Rcpp::stop("pollen kernel other than 3 specified");
	    }
	}
      //      Rprintf("dens.size() %i",dens.size());
      ds = 0;
      for (k=0;k<dens.size();k++) ds += dens[k];//Rprintf("dens %g\n",dens[k]);
      for (k=0;k<dens.size();k++) dens[k]= dens[k]/ds;//Rprintf("dens %g\n",dens[k]);
      if (dens.size()>0)
	{
	  dad= which_max(oneMultinomCall(wrap(dens)));
	  ret(i,1) = oid[dad];
	  //	  	  Rprintf("just assigned %i\n",i);

	}
      else
	{
	  ret(i,1)= -1 ;
	}
    }
  
  return ret;
}

// [[Rcpp::export]]
NumericVector FindMalesDist(IntegerVector reprows,
			   NumericMatrix individuals,
			   NumericVector mrdist //sample distance for each male
			   )

{ //8 times slower than FindMales4

  //  Timer timer;
  int i, j, k, r, dad;
  double x,y,mind;
  NumericVector mrdist2=clone(mrdist);
  NumericMatrix ret(reprows.size(),2);
  ret(_,0) = reprows;

  NumericVector dsts(individuals.nrow());

  mrdist2=pow(mrdist2,2);  //gonna compare a^2+b^2 to rdists^2 and try to avoid sqrt()
  
  for (i=0; i<reprows.size(); i++)
    {//figure out which pops are overlaid by the circle. only keep those individuals

      //timer.step("start");
      
      r = reprows[i];
      x=individuals(r,3);
      y=individuals(r,4);

      dsts = (pow(x-individuals(_,3),2)+pow(y-individuals(_,4),2));
      mind=100000000;
      dad=-1;

      for (j=0;j<dsts.size();j++)
	  if (abs(dsts[j]-mrdist2[i])<mind)
	    {
	      mind=abs(dsts[j]-mrdist[i]);
	      dad=j;
	    }
      ret(i,1)=dad+1;
    }
  return ret;
}


// [[Rcpp::export]]
NumericMatrix Inheritance(NumericMatrix omat,
			  NumericMatrix inds,
			  NumericVector lvec,
			  NumericVector pl
			  )
{

  int i, j, cnt;
  
  NumericMatrix newind(omat.nrow(),inds.ncol());
  
  NumericVector mom(inds.ncol());
  NumericVector dad(inds.ncol());
  
  for (i=0; i<omat.nrow(); i++)
    {
      mom  = inds(omat(i,0)-1,_);  //-1 is to convert from r indexinig
      dad  = inds(omat(i,1)-1,_);  //-1 is to convert from r indexinig

      cnt=9;
      for (j=0;j<pl.size();j++)
	{
	  if (pl[j]==1) //haploid is maternal only, should improve
	    {
	      newind(i,cnt)=mom[cnt];
	      cnt++;
	    }
	  else //ploidy 2, basically.  need to think about arbitrary ploidy and generalize
	    {

	      newind(i,cnt)=ifelse(runif(1)>0.5,mom[cnt],mom[cnt+1])[0];
	      newind(i,cnt+1)=ifelse(runif(1)>0.5,dad[cnt],dad[cnt+1])[0];
	      
	      cnt=cnt+2;
	    }
	}
      
    }
  return newind;
}
