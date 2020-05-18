#include <cstdint>
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
std::size_t lastIdx; //for looking at the complement
std::ofstream output("frames.txt");

//Checks for the start codon from complement sequence
std::size_t startCodonComplement(){
    lastIdx-=3; //Incrememnt idx so matches aren't repeated
    for(std::size_t firstLetterOfStartCodon;lastIdx > 2; lastIdx-=3){
        //Check if it's an ATG by looking at ATG's reverse complement CAT
        if(totalLine[lastIdx] < 'D'){
            //Then totalLine[lastIdx] must be an 'A' or 'C'
            //Must return an index that equals the 'T', since that's 'A's complement and A is first letter in ATG 
            if(totalLine[lastIdx] == 'C' && totalLine[lastIdx + 1] == 'A' 
                && totalLine[firstLetterOfStartCodon = (lastIdx + 2)] == 'T'){
                //std::cout << "first: " << totalLine[lastIdx] << totalLine[lastIdx + 1] << totalLine[lastIdx + 2] << '\n';
                return firstLetterOfStartCodon;//start codon found
            }
            //If the above if condition isn't true, then totalLine[lastIdx] == 'A'
            else if(totalLine[lastIdx - 1] == 'C' && totalLine[firstLetterOfStartCodon = (lastIdx + 1)] == 'T'){
                //std::cout << "second: " << totalLine[lastIdx - 1] << totalLine[lastIdx] << totalLine[lastIdx + 1] << '\n';
                return firstLetterOfStartCodon;//start codon found
            }
            //If neither the if or else if conditions are true, no start codon is found

        } else if(totalLine[firstLetterOfStartCodon = lastIdx] == 'T' && totalLine[lastIdx - 1] == 'A'
                 && totalLine[lastIdx - 2] == 'C') {
            //std::cout << "third: " << totalLine[lastIdx - 2] << totalLine[lastIdx - 1] << totalLine[lastIdx] << '\n';
            return firstLetterOfStartCodon;//start codon found
            //It's a T. 
        }
        //If nothing else, it's a G and no matches are possible
    }
    return std::string::npos;//Reached end of string, no matches are available
}

std::size_t stopCodonComplement(std::size_t start){
    for(std::size_t second{start - 1}, third{start - 2}; third > 2; start-=3,second-=3,third-=3){
        //output << totalLine[first] << totalLine[second] << totalLine[third]; //The sequence. The start and stop codons
        if(totalLine[start] == 'A' && ((totalLine[second] == 'T'
         && (totalLine[third] == 'C' || totalLine[third] == 'T')) 
         || (totalLine[second] == 'C' && totalLine[third] == 'T'))) {
             //stop codon found
             //std::cout << "Stop codon: " << totalLine[first] << totalLine[second] << totalLine[third] << '\n';
             return start - 3; //The last index of the sequence.
             //Sequence's length should be lastIdx - (first - 3)
        }
    }
    return 0; //No stop codons
}

std::size_t stopCodonNormal(std::size_t start){
    for(std::size_t second{start + 1}, third{start + 2}; third < totalLine.size();
        start+=3, second+=3, third+=3){
        //output << totalLine[first] << totalLine[second] << totalLine[third];
        if(totalLine[start] == 'T' && 
            ((totalLine[second] == 'A' && (totalLine[third] == 'A' || totalLine[third] == 'G')) || 
            (totalLine[second] == 'G' && totalLine[third] == 'A'))){
                return start + 3;
            }
    }
    return 0;
}

constexpr bool compareFrames(std::size_t beginOne, std::size_t endOne, uint_least16_t frameOne,
                   std::size_t beginTwo, std::size_t endTwo, uint_least16_t frameTwo){
                       //If true, frame 1 is printed to file before frame 2
    return (beginOne == beginTwo ?  
    (endOne == endTwo ? frameOne < frameTwo : endOne < endTwo) : beginOne < beginTwo);
}

void loadStringIntoFile(std::size_t start, std::size_t end){
    //for(;end > start; --end) output << totalLine[end];
    for(;start < end; ++start) output << totalLine[start];
    output << '\n';
}
//taatgtgtaaaattaattttagtagtgctatccccat
void completeMethod(){

    //The initial values of complement are in the 29000s, but you also have to change it to the low-ish equivalent, where the idx is targeting the beginning letters of a string. MAKE SURE YOU DON'T CONFUSE THE TWO
    std::size_t normal{totalLine.find(line = "ATG", 0)}, complement{startCodonComplement()};
    for(std::size_t resultStopCodonNormal{stopCodonNormal(normal)}, resultStopCodonComplement{stopCodonComplement(complement)}; normal != std::string::npos || complement != std::string::npos;){
        //const bool smaller(compareFrames(normal, resultStopCodonNormal = stopCodonNormal(normal), normal % 3, 
        //complement, resultStopCodonComplement = stopCodonComplement(), complement % 3)); 
        if constexpr(false){
            std::cout << "breakpoint\n";
        }

        if(compareFrames(normal, resultStopCodonNormal, normal % 3,totalLine.size() - complement, resultStopCodonComplement, complement % 3)){
            output << normal << '\t' << resultStopCodonNormal << '\t' << (normal % 3) << '\t';
            for(; normal < resultStopCodonNormal; ++normal) output << totalLine[normal];
            output << '\n';
            resultStopCodonNormal = stopCodonNormal(normal = totalLine.find("ATG", normal + 1));
        } else {
            output << (totalLine.size() - complement) << '\t' << (totalLine.size() - resultStopCodonComplement) << '\t' << ((complement % 3) + 4) << '\t';

            for(std::size_t idx{complement}; idx > resultStopCodonComplement; --idx){
                output << totalLine[idx];
            }
            output << '\n';
            resultStopCodonComplement = stopCodonComplement(complement = startCodonComplement());
        }
    }
}


int main(){
    
    //skip first line
    for(getline(theFile, line);getline(theFile, line); ){
        totalLine+=line;
    }
    lastIdx = totalLine.size(); //Tested with 4000 and 7500

    //Must change lastIdx's value when match is found so you don't keep reporting the same match
    /*
    std::cout << "Begin startCodonComplement\n" << startCodonComplement()
    <<  ' ' << startCodonComplement()
     <<  ' ' << startCodonComplement()
     << '\n';
     */
     completeMethod();
/*
     std::cout << "Begin stopCodonComplement\n";
     lastIdx = 600;
     stopCodonComplement();
*/
    //Print statements for debugging
    //std::cout << totalLine[106] << totalLine[106 + 1] << totalLine[106 + 2] << '\n';
    //std::cout << totalLine[133] << totalLine[133 + 1] << totalLine[133 + 2] << '\n';
    if constexpr(false){
    for(std::size_t start = totalLine.find(line = "ATG", 0);start != std::string::npos; start = totalLine.find(line,start+1)){
        //totalLine[start] == 'A' from "ATG"

        std::size_t first(start);
        /*
        for(std::size_t second(first + 1),third(first + 2); third < totalLine.size();first+=3, second+=3, third+=3){
            output << totalLine[first] << totalLine[second] << totalLine[third]; //The sequence. The start and stop codons and everything in between.
            //These ifs look for stop codons
            if(totalLine[first] == 'T' && 
            ((totalLine[second] == 'A' && (totalLine[third] == 'A' || totalLine[third] == 'G')) || 
            (totalLine[second] == 'G' && totalLine[third] == 'A'))){
                break;
            }
        }
        */
        //The char at totalLine[first],totalLine[first + 1],totalLine[first + 2] form a stop codon
        output << '\t' << start << '\t' << stopCodonNormal(first) << '\t' << (start % 3) << '\n';
        //Tab + start + tab + end + tab + frame number + new line
        //std::cout << start << ' ';
    }
    }
}
//clang++ -Wall -std=c++1z -stdlib=libc++ -g main.cpp -o main && ./main
//AAAATGCTTAAACCATTGCCC