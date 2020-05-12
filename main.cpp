#include <cstdint>
#include <iostream>
#include <fstream> //open files
#include <vector>

struct Frame{
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
    line = "ATG";
    std::ofstream output("frames.txt");
    for(std::size_t pos = totalLine.find(line, 0);pos != std::string::npos; pos = totalLine.find(line,pos+1)){
        std::size_t third(pos + 2);
        for(std::size_t first(pos), second(pos + 1); third < totalLine.size(); ++third){
            //These ifs look for stop codons
            if(totalLine[first++] == 'A'){
                if(totalLine[second] == 'T'){
                    if(totalLine[third] == 'T' || totalLine[third] == 'C'){
                        break;
                    }
                }
                else if(totalLine[second] == 'C'){
                    if(totalLine[third] == 'T'){
                        break;
                    }
                }
            }
            ++second;
        }
        output << pos << '\t' << third << '\t' << (third - pos) << '\n';
        //Here, find open reading frames. No need to use vectors
    }



/*
    for(const auto& i : positions) std::cout << i << ' ';
    std::cout << '\n';
*/
    //std::cout << totalLine[17680] << totalLine[17680 + 1] << totalLine[17680 + 2] << '\n';
    

}