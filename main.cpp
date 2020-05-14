#include <cstdint>
#include <iostream> //std::cout for print statements to debug stuff
#include <fstream> //open files

struct Frame{
    //Unused and probably unnecessary Frame object
    const std::size_t start, end, frame;
    constexpr Frame() : start{0}, end{0}, frame{0}{}
    constexpr Frame(std::size_t s, std::size_t e, std::size_t f) : start{s}, end{e}, frame{f}{}
};

int main(){
    std::ifstream theFile("sars_cov2.fasta");
    std::string line, totalLine;
    //skip first line
    for(getline(theFile, line);getline(theFile, line); ){
        totalLine+=line;
    }
    //Print statements for debugging
    //std::cout << totalLine[106] << totalLine[106 + 1] << totalLine[106 + 2] << '\n';
    //std::cout << totalLine[133] << totalLine[133 + 1] << totalLine[133 + 2] << '\n';
    std::ofstream output("frames.txt");
    for(std::size_t start = totalLine.find(line = "ATG", 0);start != std::string::npos; start = totalLine.find(line,start+1)){
        //totalLine[start] == 'A' from "ATG"
        std::size_t first(start);
        for(std::size_t second(start + 1),third(start + 2); third < totalLine.size();first+=3, second+=3, third+=3){

            output << totalLine[first] << totalLine[second] << totalLine[third]; //The sequence. The start and stop codons and everything in between.
            //These ifs look for stop codons
            if(totalLine[first] == 'T' && 
            ((totalLine[second] == 'A' && (totalLine[third] == 'A' || totalLine[third] == 'G')) || 
            (totalLine[second] == 'G' && totalLine[third] == 'A'))){
                break;
            }
            
        }
        //The char at totalLine[first],totalLine[first + 1],totalLine[first + 2] form a stop codon
        output << '\t' << start << '\t' << (first + 3) << '\t' << (start % 3) << '\n';
        //Tab + start + tab + end + tab + frame number + new line
        //std::cout << start << ' ';
    }
}
//clang++ -Wall -std=c++1z -stdlib=libc++ -g main.cpp -o main && ./main