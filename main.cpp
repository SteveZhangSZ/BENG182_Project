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
    for(std::size_t pos = totalLine.find(line = "ATG", 0);pos != std::string::npos; pos = totalLine.find(line,pos+1)){
        std::size_t first(pos);
        for(std::size_t second(pos + 1),third(pos + 2); third < totalLine.size();first+=3, second+=3, third+=3){
            //These ifs look for stop codons
            if(totalLine[first] == 'T' && 
            ((totalLine[second] == 'A' && (totalLine[third] == 'A' || totalLine[third] == 'G')) || 
            (totalLine[second] == 'G' && totalLine[third] == 'A'))){
                break;
            }
            output << totalLine[first] << totalLine[second] << totalLine[third]; //The sequence
        }
        output << '\t' << pos << '\t' << first << '\t' << (pos % 3) << '\n';
        //Tab + start + tab + end + tab + frame number + new line
    }
}
//clang++ -Wall -std=c++1z -stdlib=libc++ -g main.cpp -o main && ./main