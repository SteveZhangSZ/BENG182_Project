#include <cstdint> //for uint_least16_t
#include <iostream> //std::cout for print statements to debug stuff
#include <fstream> //open files

#define UNDERFLOWPREVENTNUM 18446744073709551612UL //Prevent underflowing to big num when working with unsigned long
struct Frame{
    //Unused and probably unnecessary Frame object
    const std::size_t start, end, frame;
    constexpr Frame() : start{0}, end{0}, frame{0}{}
    constexpr Frame(std::size_t s, std::size_t e, std::size_t f) : start{s}, end{e}, frame{f}{}
};
std::ifstream theFile("sars_cov2.fasta");
std::string line, totalLine;
std::size_t lastIdx; //for looking at the complement only
std::ofstream output("frames.txt");

//Looks for the start codon locations in complement sequence
std::size_t startCodonComplement(){
    lastIdx-=3; //Increment idx so matches aren't repeated
    for(std::size_t firstLetterOfStartCodon;lastIdx < UNDERFLOWPREVENTNUM; lastIdx-=3){ 
        //Check if it's an ATG by looking at ATG's reverse complement CAT
        //Must return an index that equals the 'T', since that's 'A's complement and A is first letter in ATG 
        if(totalLine[lastIdx] == 'C' && totalLine[lastIdx + 1] == 'A' 
            && totalLine[firstLetterOfStartCodon = (lastIdx + 2)] == 'T'){
            return firstLetterOfStartCodon;//start codon found
        }
        else if(totalLine[lastIdx] == 'A' &&
            totalLine[lastIdx - 1] == 'C' && totalLine[firstLetterOfStartCodon = (lastIdx + 1)] == 'T'){
            return firstLetterOfStartCodon;//start codon found
        }
        else if(totalLine[firstLetterOfStartCodon = lastIdx] == 'T' && totalLine[lastIdx - 1] == 'A'
                 && totalLine[lastIdx - 2] == 'C') {
            return firstLetterOfStartCodon;//start codon found
        }
        //If nothing else, it's a G and no matches are possible
    }
    return std::string::npos;//Reached end of string, no matches are available
}
//Looks for the stop codon locations in complement sequence
std::size_t stopCodonComplement(std::size_t start){
    for(std::size_t second{start - 1}, third{start - 2}; third < UNDERFLOWPREVENTNUM; start-=3,second-=3,third-=3){
        if(totalLine[start] == 'A' && ((totalLine[second] == 'T'
        && (totalLine[third] == 'C' || totalLine[third] == 'T')) 
        || (totalLine[second] == 'C' && totalLine[third] == 'T'))){
            return start - 3; //The last index of the sequence. Stop codon found
            //Sequence's length should be lastIdx - (first - 3)
        }
    }
    return 0; //No stop codons
}
//Looks for the stop codon locations in normal/forward sequence
std::size_t stopCodonNormal(std::size_t start){
    for(std::size_t second{start + 1}, third{start + 2}; third < totalLine.size();
        start+=3, second+=3, third+=3){
        if(totalLine[start] == 'T' && 
            ((totalLine[second] == 'A' && (totalLine[third] == 'A' || totalLine[third] == 'G')) || 
            (totalLine[second] == 'G' && totalLine[third] == 'A'))){
                return start + 3;
            }
    }
    return totalLine.size();
}
//Compares the frames to decide which to add to the file first
constexpr bool compareFrames(std::size_t beginOne, std::size_t endOne, uint_least16_t frameOne,
                   std::size_t beginTwo, std::size_t endTwo, uint_least16_t frameTwo){
    return (beginOne == beginTwo ?  //If method returns true, frame 1 is printed to file before frame 2
    (endOne == endTwo ? frameOne < frameTwo : endOne < endTwo) : beginOne < beginTwo);
}
//The method to write into the file
void completeMethod(){
    std::size_t normal{totalLine.find(line = "ATG", 0)}, complement{startCodonComplement()};
    for(std::size_t resultStopCodonNormal{stopCodonNormal(normal)}, 
    resultStopCodonComplement{stopCodonComplement(complement)}; 
    normal != std::string::npos || complement != std::string::npos;){
        if(compareFrames(normal, resultStopCodonNormal, normal % 3,totalLine.size() - complement, 
        totalLine.size() - resultStopCodonComplement, complement % 3)){
            if(resultStopCodonNormal - normal >= 90){ //Checks if there are 90 nucleotides or more
                output << normal << '\t' << resultStopCodonNormal << '\t' << (normal % 3) << '\t';
                for(; normal < resultStopCodonNormal; ++normal) output << totalLine[normal];
                output << '\n';
            }
            resultStopCodonNormal = stopCodonNormal(normal = totalLine.find("ATG", normal + 1));
        } else {
            if(complement - resultStopCodonComplement >= 90){
                output << (totalLine.size() - complement) << '\t' << (totalLine.size() - resultStopCodonComplement) << '\t' << ((complement % 3) + 4) << '\t';
                for(std::size_t idx{complement}; idx > resultStopCodonComplement; --idx){
                    output << 
                    (totalLine[idx] < 'D' ? (totalLine[idx] == 'A' ? 'T' : 'G') : (totalLine[idx] == 'T' ? 'A' : 'C'));
                    //Puts the complement into the file
                }
                output << '\n';
            }
            resultStopCodonComplement = stopCodonComplement(complement = startCodonComplement());
        }
    }
}
int main(){
    //skip first line
    for(getline(theFile, line);getline(theFile, line); ){
        totalLine+=line;
    }
    lastIdx = totalLine.size();
    completeMethod();
}
//clang++ -Wall -std=c++1z -stdlib=libc++ -g main.cpp -o main && ./main