#include <cstdint> //for uint_least16_t
#include <iostream> //std::cout for print statements to debug stuff
#include <fstream> //open files
struct Frame{
    //Unused and probably unnecessary Frame object
    const std::size_t start, end, frame;
    constexpr Frame() : start{0}, end{0}, frame{0}{}
    constexpr Frame(std::size_t s, std::size_t e, std::size_t f) : start{s}, end{e}, frame{f}{}
};
std::ifstream theFile("sars_cov2.fasta");
std::string line, totalLine;
std::size_t startCodonNormal(std::size_t start){
    for(/*std::size_t firstLetterOfStartCodon*/; start < totalLine.size(); start+=3){
    /*
        if(start == 106){
            std::cout << totalLine[start] << totalLine[start + 1] << totalLine[start + 2] << '\n';
        }
    */
        if(totalLine[start] == 'A' && totalLine[start + 1] == 'T'
        && totalLine[start + 2] == 'G'){
            return start;//start codon found
        }
        /*
        if(totalLine[start] == 'G' && totalLine[start - 1] == 'T'
        && totalLine[firstLetterOfStartCodon = start - 2] == 'A'){
            return firstLetterOfStartCodon;//start codon found
        }
        if(totalLine[start] == 'T' && totalLine[firstLetterOfStartCodon = start - 1] == 'A'
        && totalLine[start + 1] == 'G'){
            return firstLetterOfStartCodon;//start codon found
        }
        */
    }
    return std::string::npos;//Reached end of string, no matches are available
}

std::size_t lastIdx{0}; /*for looking at the complement only*/

/*
1. Initialize a start codon value for a frame

Loop:
2. Find a stop codon that's reachable by the same frame 
and greater than the current start codon value
3. Change the same frame's start codon index value to be an 
appropriate one that's greater than the stop codon's index
*/
std::ofstream output("frames.txt");

//Start == 2, 3, or 4, different value for each frame. Start idx also equals the first letter of the codon
std::size_t findFirstStopCodonComplement(std::size_t start){
    for(; start < totalLine.size(); start+=3){
        if(totalLine[start] == 'A' && ((totalLine[start - 1] == 'T' && (totalLine[start - 2] == 'C'
        || totalLine[start - 2] == 'T') ) || (totalLine[start - 1] == 'C' && totalLine[start - 2] == 'T')) ){
            return start;
        }
    }
    return std::string::npos;
}

/* How findTwoStopCodons should work:

Finds the first and second stop codon on a frame
Then second and third on the same frame
Then third and fourth etc.
 
The very first value of its parameter, start, will be the index of the
Very first stop codon of that frame. 
 
It finds start codons, only keeping the maximum index, until it reaches another stop codon.
 
The returned value, dubbed RV, will be the “another stop codon” and be the 
start parameter upon the next call
*/
std::pair<std::size_t,std::size_t> findStartCodonBetweenStop(std::size_t start){//For the reverse complement
    //Find the indices of two subsequent stop codons on the same frame
    for(std::size_t theStart{0}; start < totalLine.size(); start+=3){
        if(totalLine[start] == 'T' && totalLine[start - 1] == 'A' &&
        totalLine[start - 2] == 'C'){ //A start codon!
            theStart = start;
        }
        else if(totalLine[start] == 'A' && ((totalLine[start - 1] == 'T' && (totalLine[start - 2] == 'C'
        || totalLine[start - 2] == 'T') ) || (totalLine[start - 1] == 'C' && totalLine[start - 2] == 'T')) ){
            return {std::move(theStart), start};
            //reached a stop codon. The start codon should be between start, 
            //the position of the first stop codon, and this start + 2 value, the second codon's position.
            //There's a possibility no start codon exists between 2 stop codons
        }
    }
    return {std::string::npos,std::string::npos};
}
//Look for stop codons in the same frame and whose index exceeds "start"


//Looks for the stop codon locations in normal/forward sequence
std::size_t otherStopCodonNormal(std::size_t start){
    //std::cout << start << ' ';
    for(std::size_t second{start + 1}, third{start + 2}; third < totalLine.size();
        start+=3, second+=3, third+=3){
        if(totalLine[start] == 'T' && 
            ((totalLine[second] == 'A' && (totalLine[third] == 'A' || totalLine[third] == 'G')) || 
            (totalLine[second] == 'G' && totalLine[third] == 'A'))){
                return start; //DO NOT include the stop codon when writing to file
            }
    }
    return totalLine.size();
}
void otherCompleteMethod(){
    std::size_t startCodonArray[]{0,0,0}; //Has all three forward frame's start index

    //Has all 6 frames', forward and reverse, stop codon index
    std::size_t stopCodonArray[]{otherStopCodonNormal(startCodonArray[0] = startCodonNormal(0)),
    otherStopCodonNormal(startCodonArray[1] = startCodonNormal(1)),
    otherStopCodonNormal(startCodonArray[2] = startCodonNormal(2)),
    findFirstStopCodonComplement(2), findFirstStopCodonComplement(3),findFirstStopCodonComplement(4)};

    //The .first has the start codon's index, which is greater than findFirstStopCodonComplement()'s
    //The .second has the stop codon's index that's further to the right of .first stop codon
    std::pair<std::size_t,std::size_t> startAndStopArray[]{findStartCodonBetweenStop(stopCodonArray[3] + 3)
    ,findStartCodonBetweenStop(stopCodonArray[4] + 3),findStartCodonBetweenStop(stopCodonArray[5] + 3)};

    for(const uint_least8_t i : {0,1,2}){
        while(!startAndStopArray[i].first){
            stopCodonArray[i + 3] = startAndStopArray[i].second;
            startAndStopArray[i] = findStartCodonBetweenStop(startAndStopArray[i].second + 3);
        } 
    }
    uint_least8_t normal(startCodonArray[0] < startCodonArray[1] ? 
    (startCodonArray[0] < startCodonArray[2] ? 0 : 2) :
    (startCodonArray[1] < startCodonArray[2] ? 1 : 2)), 
    complement(startAndStopArray[0].first < startAndStopArray[1].first ? 
    (startAndStopArray[0].first < startAndStopArray[2].first ? 0 : 2) :
    (startAndStopArray[1].first < startAndStopArray[2].first ? 1 : 2)); //The values are indices of arrays

    //Now normal and complement have the indices of the smallest forward and reverse starting values.

    uint_least8_t compPlusThree(complement + 3);

    //Must make sure all start codon values don't equal std::string::npos
    while(startCodonArray[0] != std::string::npos || startCodonArray[1] != std::string::npos ||
    startCodonArray[2] != std::string::npos || startAndStopArray[0].first != std::string::npos ||
    startAndStopArray[2].first != std::string::npos || startAndStopArray[1].first != std::string::npos){
        //Note the <=. Comparing start codon indices here
        if(startCodonArray[normal] <= startAndStopArray[complement].first){ 
        
            if(stopCodonArray[normal] - startCodonArray[normal] >= 90){
                output << startCodonArray[normal] << '\t' << (stopCodonArray[normal] + 1) << "\t+\t";
                for(std::size_t idx{startCodonArray[normal]}; idx < stopCodonArray[normal]; ++idx){
                    output << totalLine[idx];
                }
                output << '\n';
                //record info into the file
            }
            stopCodonArray[normal] = otherStopCodonNormal(startCodonArray[normal] = 
            startCodonNormal(stopCodonArray[normal] + 3));
            normal = (startCodonArray[0] < startCodonArray[1] ? 
            (startCodonArray[0] < startCodonArray[2] ? 0 : 2) :
            (startCodonArray[1] < startCodonArray[2] ? 1 : 2));
        } else {
            if(startAndStopArray[complement].first - stopCodonArray[compPlusThree] >= 90){
                output << startAndStopArray[complement].first << '\t' << (stopCodonArray[compPlusThree] - 1) << "\t-\t";
                for(std::size_t idx{startAndStopArray[complement].first}; idx > stopCodonArray[compPlusThree]; --idx){
                    output << 
                    (totalLine[idx] < 'D' ? (totalLine[idx] == 'A' ? 'T' : 'G') : (totalLine[idx] == 'T' ? 'A' : 'C'));
                    //Puts the complement into the file
                }
                output << '\n';
                //record info into the file
            }
            do{
                startAndStopArray[complement] = findStartCodonBetweenStop((stopCodonArray[compPlusThree]
                 = startAndStopArray[complement].second) + 3);
            } while(!startAndStopArray[complement].first);
            compPlusThree = (complement = (startAndStopArray[0].first < startAndStopArray[1].first ?
            (startAndStopArray[0].first < startAndStopArray[2].first ? 0 : 2) :
            (startAndStopArray[1].first < startAndStopArray[2].first ? 1 : 2))) + 3;
        }
    }
}
int main(){
    //skip first line
    for(getline(theFile, line);getline(theFile, line); ){
        totalLine+=line;
    }
    lastIdx = totalLine.size();
    otherCompleteMethod();
}
//clang++ -Wall -std=c++1z -stdlib=libc++ -g main.cpp -o main && ./main