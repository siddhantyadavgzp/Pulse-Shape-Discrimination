/*-- 
The following code does numerical integration on the wave data to 
product time constant and charge before speration line and after seperation line. 
--*/
/*-- Input Format --*/
/*-- 
    No_of_pulses
    pulse[1]
    pulse[2]
    .
    .
    .
    pulse[no]
--*/
/*--  Output Format --*/
/*-- 

    Sl_no,q1,q2,_constant,Classification

--*/
     

#include "bits/stdc++.h"

using namespace std;




bool fo=0;

void solve(){
    int n;
    cin>>n;
    for ( int i=0;i<n;i++){
        vii va;
        int mipos=0;
        double  tot=0;
        ll mx=900;
        for ( int j=0;j<191;j++){
            double xx;
            cin>>xx;
            va.pb( xx);
            if ( va[j] < va[mipos])
                mipos=j;
        }
        assert( sz( va)==191);
        for ( int j=0;j<sz(va);j++)
            va[j]=mx-va[j];
        int mark=mipos+17;
        for ( int j=mipos;j<191;j++){
            if ( va[j]*2<va[mipos]){
                mark=j;
                break;
            }
        }

        cout<<mark-mipos<<',';
        mark=57;
        for ( int j=0;j<mark;j++){
            tot+=va[j];
        }
        cout<<tot<<',';
        tot=0;
        for ( int j=mark;j<191;j++){
            tot+=va[j];
        }
        cout<<tot<<","<<fo<<"\n";
    }
}

int main(){

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    #ifndef ONLINE_JUDGE
        if ( fo){
            // Input file for gamma
            freopen("FGamma.txt","r",stdin);
        }
        else{
            // Input file for neutron
            freopen("FNeutron.txt","r",stdin);
        }
        freopen("OUTPUT.txt","w",stdout);
    #endif

    solve()
}