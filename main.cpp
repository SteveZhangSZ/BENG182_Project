#include <cstdint>
#include <iostream>
#include <fstream> //open files

struct Frame{
    //Unused Frame object
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
    line = "ATG"; //the start codon
    //Print statements for debugging
    //std::cout << totalLine[106] << totalLine[106 + 1] << totalLine[106 + 2] << '\n';
    //std::cout << totalLine[133] << totalLine[133 + 1] << totalLine[133 + 2] << '\n';
    std::ofstream output("frames.txt");
    for(std::size_t pos = totalLine.find(line, 0);pos != std::string::npos; pos = totalLine.find(line,pos+1)){
        std::size_t first(pos);
        for(std::size_t second(pos + 1),third(pos + 2); third < totalLine.size();third = (second = (first+=3) + 1) + 1){
            //These ifs look for stop codons
            if(totalLine[first] == 'T'){
                if(totalLine[second] == 'A' && (totalLine[third] == 'A' || totalLine[third] == 'G')){
                    break;
                }
                else if(totalLine[second] == 'G' && totalLine[third] == 'A'){
                    break;
                }
            }
        }
        output << pos << '\t' << first << '\t' << (first - pos) << '\n';
        //Here, find open reading frames. No need to use vectors
    }
}