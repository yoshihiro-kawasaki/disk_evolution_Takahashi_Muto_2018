#include "utils.hpp"

#include <iostream>
#include <fstream>

namespace Utils {


double LinearInterPolation(double x, double x1, double x2, double y1, double y2)
{
    if (x1 > x2) {
        std::cerr << "Utility::LinearInterPolation : error (x1 > x2) " << std::endl;
        std::exit(1);
    }

    return (y1 + ((( y2 - y1 ) * ( x - x1 )) / ( x2 - x1 )));
}

std::vector<std::string> Split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    std::string item;
    for (char ch: s) {
        if (ch == delim) {
            if (!item.empty())
                elems.push_back(item);
            item.clear();
        }
        else {
            item += ch;
        }
    }
    if (!item.empty())
        elems.push_back(item);
    return elems;
}

void FileCopy(const std::string filename_in, const std::string filename_out)
{
    std::ifstream file_in(filename_in, std::ios::in);
    if(!file_in) {
        std::cerr << "##ERROR : Failed to open " << filename_in << std::endl; exit(1); 
    }
    std::ofstream file_out(filename_out, std::ios::trunc | std::ios::out);
    if(!file_out) {
        std::cerr << "##ERROR : Failed to open " << filename_out << std::endl; exit(1); 
    }

    std::string str;

    while (std::getline(file_in, str)) {
        file_out << str << std::endl;
    }

    file_in.close();
    file_out.close();
    
    return;
}

bool ConvertStrToBool(std::string strbool)
{
    if (strbool == "true" || strbool == "True" || strbool == "TRUE") {
        return true;
    } else if (strbool == "false" || strbool == "False" || strbool == "FALSE") {
        return false;
    } else {
        std::cerr << "## Error : Utils::ConvertStrToBool, str = " << strbool << std::endl;
        std::cerr << "str must be 'true' or 'false'." << std::endl;
        std::exit(1);
    }
}

// class InputConfigure

InputConfigure::InputConfigure()
{

}

InputConfigure::InputConfigure(std::string file_name)
{
    std::ifstream file(file_name, std::ios::in);
    if (!file) {
        std::cout << "error open file : " << file_name << std::endl;
        exit(1);
    }

    input_filename = file_name;

    std::string str;
    std::vector<std::string> line;

    while (std::getline(file, str)) {
        std::cout << str[0] << std::endl;
        if (str[0] == '#') continue;
        
        line = Split(str, ' ');
        dict[line[0]] = line[2];
    } 

    file.close(); 
}

InputConfigure::~InputConfigure()
{

}

bool InputConfigure::Contains(const std::string &key)
{
    return dict.find(key) != dict.end();
}

void InputConfigure::Insert(std::string key, std::string item)
{
    dict[key] = item;
}

void InputConfigure::ReadFile(std::string file_name)
{
    std::ifstream file(file_name, std::ios::in);
    if (!file) {
        std::cout << "error open file : " << file_name << std::endl;
        exit(1);
    }

    input_filename = file_name;

    std::string str;
    std::vector<std::string> line;

    while (std::getline(file, str)) {
        // std::cout << str << std::endl;
        if (str[0] == '#' || str.empty()) continue;
        
        line = Split(str, ' ');
        dict[line[0]] = line[2];
    } 

    file.close();

    return;
}

size_t InputConfigure::Size()
{
    return dict.size();
}

std::string InputConfigure::GetString(std::string key)
{
    if (Contains(key)) {
        return dict[key];
    } else {
        std::cout << "## No item : " << key << std::endl;
        return "";
    }
}

double InputConfigure::GetDouble(std::string key)
{
    if (Contains(key)) {
        return std::stod(dict[key]);
    } else {
        std::cout << "## No item : " << key << std::endl;
        return 0.0;
    }
}


int InputConfigure::GetInt(std::string key)
{
    if (Contains(key)) {
        return std::stoi(dict[key]);
    } else {
        std::cout << "## No item : " << key << std::endl;
        return 0;
    }
}


bool InputConfigure::GetBool(std::string key)
{
    if (Contains(key)) {
        return Utils::ConvertStrToBool(dict[key]);
    } else {
        std::cout << "## No item : " << key << std::endl;
        return 0;
    }
}

std::string InputConfigure::GetInputFileName()
{
    return input_filename;
}

} // end namespace Utils