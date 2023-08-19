#include <Rcpp.h>
#include <vector>

using namespace Rcpp;





// Other functions used
//###################################################################################################################

// Simple function to display a 1-d vector
// template <typename T>
// void display(std::vector<T> const &vector){
//   
//   for (int i = 0; i < vector.size(); i++)
//     std::cout <<std::fixed << std::setprecision(2) << std::setw(5) << std::setfill(' ') << vector[i] << ' ';
//   
//   std::cout<<"\n";
// }




template <typename T>
std::vector<size_t> sort_indexes(std::vector<T> const &v)
{
  
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  
  return idx;
}





template <typename T>
std::vector<double> calculate_ranks(const std::vector<T>& v, std::vector<double>& countOfElements) {
  
    std::vector<T> sortedCopy = v; // Make a copy of the original vector
    std::sort(sortedCopy.begin(), sortedCopy.end()); // Sort the copy to obtain ranks
    
    std::vector<double> ranks(v.size());
    std::unordered_map<T, double> sumRanks;
    std::unordered_map<T, size_t> countElements;
    
    for (size_t i = 0; i < sortedCopy.size(); ++i) {
      sumRanks[sortedCopy[i]] += i;
      countElements[sortedCopy[i]]++;
    }
    // // iterating over all value of sumRanks
    //   typename unordered_map<T, double>::iterator itr;
    //   cout << "\nsumRanks : \n";
    //   for (itr = sumRanks.begin(); 
    //        itr != sumRanks.end(); itr++) 
    //   {
    //     // itr works as a pointer to 
    //     // pair<string, double> type 
    //     // itr->first stores the key part and
    //     // itr->second stores the value part
    //     cout << itr->first << "  " << 
    //             itr->second << endl;
    //   }
    
    
    
    //   // iterating over all value of countElements
    //   typename unordered_map<T, size_t>::iterator itr2;
    //   cout << "\ncountElements : \n";
    //   for (itr2 = countElements.begin(); 
    //        itr2 != countElements.end(); itr2++) 
    //   {
    //     // itr works as a pointer to 
    //     // pair<string, double> type 
    //     // itr->first stores the key part and
    //     // itr->second stores the value part
    //     cout << itr2->first << "  " << 
    //             itr2->second << endl;
    //   }
    
    
    double currentRank = 1.0;
    for (size_t i = 0; i < v.size(); ++i) {
      double avgRank = sumRanks[v[i]] / countElements[v[i]];
      ranks[i] = avgRank + 1.0; // Add 1 to adjust for 1-based indexing (optional)
      }
    
    //display(ranks);
    
    
    
    for (const auto& pair : countElements) {
      if(pair.second!=1) 
        countOfElements.push_back(pair.second);
    }
    
    // std::cout << "Values in countElements:\n";
    // for (const auto& value : countElements) {
    //     std::cout << value << " ";
    // }
    // std::cout << std::endl;
    
    
    
    return ranks;
}







void mannwhitneyutest( 
    Rcpp::NumericVector const &V1,
    Rcpp::NumericVector const &V2,
    double& bothtails,
    double& righttail,
    double& lefttail
)
{
  int V1_size = V1.size();
  //copy V1 and V2 into X
  std::vector <double> X(V1.size());
  std::copy(V1.begin(), V1.end(), X.begin());
  X.insert( X.end(),V2.begin(),V2.end());
  
  
  
  int n = X.size();
  
  
  std::vector <double> ranks (X.size());
  std::vector <double> countOfElements;   // This vector will store the number of ties for the K-th rank. (size(countOfElements) <  size(ranks))
  // we use countOfElements to account for tie correction when calculating sigmaU
  
  
  // fast version of getranks  
  ranks = calculate_ranks(X, countOfElements);
  
  
  // std::cout<<"ranks : ";
  // display(ranks);
  
  // std::cout<<"countOfElements : ";
  // display(countOfElements);
  
  
  
  //ranked sums
  double T0 = std::accumulate(ranks.begin(), ranks.begin()+V1_size, 0);
  double T1 = std::accumulate(ranks.begin()+V1_size, ranks.end(), 0);
  
  
  
  
  int n0 = V1_size, n1 = X.size()-V1_size;
  
  //now i can restore X to the original size
  //X.erase(X.end()-n1, X.end());
  
  // std::cout<<"n0: "<<n0<<"\n";
  // std::cout<<"n1: "<<n1<<"\n";
  // std::cout<<"T0: "<<T0<<"\n";
  // std::cout<<"T1: "<<T1<<"\n\n";
  
  //Calculating U values:
  
  double U0, U1;
  
  U0 = n0*n1 + (n0*(n0+1)/(2))-T0;    //Per wikipedia U0 = T0-((n0*(n0+1)/2))
  U1 = n0*n1 + (n1*(n1+1)/(2))-T1;    //Per wikipedia U1 = T1-((n1*(n1+1)/2))
  //n0*n1 = U0+U1 -> sempre detto da wikipedia
  
  double U; // U= min(U0, U1)
  U0<U1 ? U = U0 : U = U1;
  // std::cout<<"U: "<<U<<"\n\n";
  
  double muU = (double(n0)*double(n1))/2;; //Expected value of U
  // std::cout<<"muU: "<<muU<<"\n";
  
  
  
  //sigmaU = (sqrt(n0)*sqrt(n1)*sqrt((n0+n1+1)))/sqrt(12);  //NOT ADJUSTED FOR TIE CORRECTION!!!
  
  
  double correctionFactor = 0;
  for (int i=0; i<countOfElements.size(); i++){
    correctionFactor += pow(countOfElements[i],3) - countOfElements[i];
  }
  
  // cout<<"n. of ties:"<<size(countOfElements)<<"\n";
  // display(countOfElements)
  
  //standard error of U
  long double sigmaU = ((sqrt(n0)*sqrt(n1)) / sqrt(12))  *   sqrt(((n0+n1+1)-(  (correctionFactor)/((n0+n1)*(n0+n1-1))   )   ));
  // std::cout<<"sigmaU: "<<sigmaU<<"\n";
  
  
  //Z value
  double z = (U-muU)/sigmaU;
  // std::cout<<"Z value: "<<z<<"\n";
  
  
  // returning the p values in the passed variables
  bothtails = (1+erf(z/sqrt(2)));  
  
  lefttail = (1+erf(z/sqrt(2)))/2;
  
  righttail = 1-lefttail;
  
}






//#############################################################################################################
// This is the actual function visible by R

// [[Rcpp::export]]
Rcpp::NumericVector fastWilcox(NumericVector V1, NumericVector V2, bool verbose = 0) {
  
  
  double bothTails = -3;
  double rightTail = -3;
  double leftTail  = -3;
  
  Rcpp::NumericVector results (3);
  
  mannwhitneyutest(V1, V2, bothTails, rightTail, leftTail);
  
  results[0] = bothTails;
  results[1] = rightTail;
  results[2] = leftTail;
  
  if(verbose){
    Rcout<<"\n";
    Rcout<<"Mann Whitney U test results: \n";
    Rcout<<"P-value bothTails: "<<bothTails<<"\n";
    Rcout<<"P-value rightTail: "<<rightTail<<"\n";
    Rcout<<"P-value leftTail : "<<leftTail<<"\n";
  }
  
  return results;
  
}










/*** R
#rcpp_MWU_Lor(c(45, 33, 35, 39, 42), c(34, 36, 41, 43, 44, 37))
*/
