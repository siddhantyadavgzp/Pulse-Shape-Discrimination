#include "bits/stdc++.h"
using namespace std;

typedef long long ll;

#define pb push_back 
#define fi first
#define se second 
#define sz(x) ((ll)x.size())
 

const double eps=1e-8;

/* 
    Structure to hold all properties realated to a pulses
    
    Contains 3 Main features 
    1) cr : used to denote charge ratio at particular energy peak
    2) peak : peak energy
    3) tp : classification

    Some other features used are
    
    1) ind: to point out the index
    2) cls: to hold the classification after the operations are performed


*/
struct st{
    int ind;
    double cr;
    int peak, tp , cls;
};


// Sorting based on peak value 
bool comp_peak( st a, st b ){    
    if ( a.peak==b.peak)
        return a.ind< b.ind;
    return a.peak < b.peak;
}   

// Sorting based on charge ratio value
bool comp_cr(st a, st b){
    return a.cr < b.cr;
}

// Function to evaluate the value of a point on a polynomial provided.
double evaluate( vector< double > &v,double x){
    double cal=1;
    double tot=0;
    for ( int i=sz(v)-1;i>=0;i--){
        tot += cal*v[i];
        cal*=x;
    }
    return tot;    
}

// Using order 4 equation for improving accuracy when count of neutorns is low
void add_cof( double a,vector< double > &v ){
    v[0]+=1;
    v[1]-=4*a;
    v[2]+=6*a*a;
    v[3]-=4*a*a*a;
    v[4]+=a*a*a*a;
}
void sub_cof( double a, vector< double > &v){
    v[0]-=1;
    v[1]+=4*a;
    v[2]-=6*a*a;
    v[3]+=4*a*a*a;
    v[4]-=a*a*a*a;
}

// Function to find the point at which the two seperate centroids yeild minima 

/*

    This function uses binary search on derative of summation [i] ( x-a_i )^4  to determine the midpoint.
    Then it is compared for both the equation and the lowest value is used.

*/
pair< double , double >  cal_min( vector< double > v, double s, double e){

    vector< double > va(sz(v)-1);
    for ( int i=0;i<4;i++){
        va[i]=v[i]*(4-i);
    }
    while ( s+eps< e){
        double mid=( s+e)/2;
        double x=evaluate( va,mid);
        if ( x>=0)
            e=mid;
        else
            s=mid;
    }



    return {s,evaluate( v,s)};
}



// To check if the dataset only contains of gamma particles and no neutrons

bool  density_check( double dcal_a, double dcal_b, vector< st > &vx,double sd){

    double m_dcal=(dcal_a+dcal_b)/2;
    int l=0,m=0, r=0;

    for ( int i=0;i<sz(vx);i++){
        if (vx[i].cr>dcal_a-sd&& vx[i].cr<dcal_a+sd){
            l++;
        }
        if (vx[i].cr>m_dcal-sd&& vx[i].cr<m_dcal+sd){
            m++;
        }
        if (vx[i].cr>dcal_b-sd&& vx[i].cr<dcal_b+sd){
            r++;
        }
    }

    if ( m<max( l,r))
        return true;

    return false;
}


// Function to sanitize data and classification

void sol( vector< st > &v){
    // Sorting based on cr value
    sort(v.begin(),v.end(), comp_cr);
    


/*
    Calculating menan and standard deviation.
    Used to ignore outliars which are 0.1 sd away from mean.
*/    

    double size_of_data=sz(v);
    double xr=0;
    for ( int i=0 ; i < sz(v) ; i++ )
        xr+=v[i].cr;
    xr/=size_of_data;
    double vr=0;
    for ( int i=0;i<sz( v);i++){
        vr+=(v[i].cr-xr)*(v[i].cr-xr);
    }
    vr=vr/size_of_data;
    double sr=sqrt(vr);
    //  sr is standard deviation and xr is the mean 

    // Removing highly arbitary data
    vector< st> vx,vy;
    for ( int i=0;i<sz(v);i++){
        if ( v[i].cr>xr-3*sr&&v[i].cr<xr+3*sr)
            vx.pb( v[i]);
        else
            vy.pb( v[i]);
    }
/*
    Creating to the power of 4 coeificients
    All of the data between range of 3 SD  are sent for classification.   

*/
    vector< double > coef(5,0),corf(5,0);

    for ( auto i: vx){
        add_cof( i.cr, coef);
    }   
    double mi=1e9,dcal_a,dcal_b;
    int sep_loc;  
    for ( int i=0;i<sz( vx)-1;i++) {
        add_cof(vx[i].cr,corf);
        sub_cof(vx[i].cr,coef);
        pair< double, double > p1=cal_min(coef,vx[i+1].cr,vx[sz(vx)-1].cr);
        pair< double, double > p2=cal_min(corf,vx[0].cr,vx[i].cr);
      //  cout<<p1.se+p2.se<<'\n';
        if ( p1.se+ p2.se < mi)
        {
            mi=p1.se+p2.se;
            sep_loc=i;
            dcal_a=p1.fi;
            dcal_b=p2.fi;
        }
    }


/* 

    Density check at midpoint
    to ensure that the data containid at centroids has more density
    Otherwise the clustering was not done at pileups rather was done on a some different curve possible gaussian


    Applying this increases accuracy by 0.1%

*/   
    if ( density_check(dcal_a,dcal_b,vx,sr*0.25)){
        for ( int i=0,j=0,k=0;i<sz(vx)||j<sz(vy);k++){
            if ( i<sz(vx) &&j<sz(vy)){
                if ( vy[i].cr<vx[sz(vx)-1].cr){
                    v[k]=vy[j++];
                    v[k].cls=0;
                }
                else if ( vy[i].cr>vx[0].cr){
                    v[k]=vy[j++];
                    v[k].cls=1;
                }
                else
                if (i<=sep_loc){
                    v[k]=vx[i++];
                    v[k].cls=0;
                }
                else{
                    v[k]=vx[i++];
                    v[k].cls=1;
                }
            }
            else
            if ( i==sz(vx)){
                v[k]=vy[j++];
                v[k].cls=1;
            }
            else if ( j==sz(vy)){
                if (i<=sep_loc){
                    v[k]=vx[i++];
                    v[k].cls=0;
                }
                else{
                    v[k]=vx[i++];
                    v[k].cls=1;
                }
            }
        }
    }
    else{
        for ( int i=0;i<sz(v);i++)
            v[i].cls=0;
    }


    sort( v.begin(),v.end(), comp_peak);

}

// Function to update the values in main array of values
void update( vector< st> &va, vector< st > & vb, int pos, int epos){
    for ( int i=pos ; i<epos ; i++ )
        va[i].cls=vb[i-pos].cls;
}

/*
    Function use case

    1) Taking input 
    2) Calling subfunctions for classification
    3) Calcualting Accuracy
    4) Printing Results

*/
void solve(){

    int n;
    cin>>n;
    
//  Taking input of the pulse data
    vector< st > va;
    for ( int i=0;i<n;i++){
        double x;
        int y,z;
        cin>>x>>y>>z;
        va.pb( {i,x, y, z,-1});
    }
    
//  Sorting based on peak values
    sort( va.begin(),va.end(), comp_peak);



//  Filtering using various analysis techniques mainly clustering.
//  Sections of energy peaks of size 50 being used to cluster out the energy


    vector< st > vb;
    int endpos=150;
    int st=0;
    for ( int i=0;i<=n;i++){
        if ( i == n || va[i].peak> endpos){    

//  Filtering function "void sol()" is called for filtering out the sepearated data. 

            sol( vb );
           
            update( va, vb, st, i );
            st=i;
            vb.clear();
            endpos += 50;
        }

        if ( i != n ) 
            vb.pb( va[i]);
    }



// Checking for accuracy
/*

    gg: True Gamma Classified as Gamma
    gn: True Neutron Classified as Gamma
    ng: True Gamma Classified as Neutron
    gg: True Neutron Classified as Neutron

*/
    int gg=0,gn=0,ng=0,nn=0;
    for ( int i=0;i<n;i++){
        if ( va[i].tp==0&&va[i].cls==0)
            gg++;
        if ( va[i].tp==0&&va[i].cls==1)
            gn++;
        if ( va[i].tp==1&&va[i].cls==0)
            ng++;
        if ( va[i].tp==1&&va[i].cls==1)
            nn++;
    }


// Calcualting accuracy
    cout<<"Gamma Classified as Gamma : "<<gg<<'\n';
    cout<<"Gamma Classified as Neutorn : "<<gn<<'\n';
    cout<<"Neutron Classified as Gamma : "<<ng<<'\n';
    cout<<"Neutron Classified as Neutron :"<<nn<<'\n';

    double num_of_sampels=gg+gn+ng+nn;
    double false_g=ng;
    double false_n=gn;
    cout<<fixed;
    cout<<setprecision(3);
    cout<<"Accuracy : " << (num_of_sampels-false_n-false_g)/num_of_sampels*100<<"%"<<'\n';




}


// Main function
int main(){

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
// File input
    freopen("INPUT.txt","r",stdin);
    freopen("OUTPUT.txt","w",stdout);

    solve();

// Printing time elapsed 
    cout <<"Time take for classification computation : "<< clock() / double(CLOCKS_PER_SEC)<<" seconds" << endl;
}
