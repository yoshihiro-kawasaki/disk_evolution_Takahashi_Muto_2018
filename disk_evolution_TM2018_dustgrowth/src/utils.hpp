#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace Utils 
{

double LinearInterPolation(double x, double x1, double x2, double y1, double y2);

/* 文字列sをある文字delimで区切って分割する */
std::vector<std::string> Split(const std::string &s, char delim);

/* ファイルをコピー */
void FileCopy(const std::string filename_in, const std::string filename_out);

/* 文字列True/Falseをboolに変換する */
bool ConvertStrToBool(std::string strbool);

// template <class T>
/* 辞書クラス */
class Dictionary
{
public:

    /* keyが辞書内に含まれているかどうかを判定 */
    bool Contains(const std::string &key) const { return dict.find(key) != dict.end(); };

    /* 辞書に新たな要素を加える */
    void Insert(std::string key, double item) { dict[key] = item; };

    /* 辞書のサイズを獲得 */
    size_t Size() { return dict.size(); }
    
    double Get_double(std::string key) 
    { 
        if (Contains(key)) {
            return  dict[key];
        } else {
            return 0.00;
        }
    }

private:
    std::map<std::string, double> dict;
    // bool valid;
};

class InputConfigure
{
public:

    InputConfigure();
    InputConfigure(std::string file_name);
    ~InputConfigure();

    bool Contains(const std::string &key);
    void Insert(std::string key, std::string item);
    void ReadFile(std::string file_name);
    size_t Size();

    std::string Get(std::string key);

    std::string GetInputFileName();

private:

    std::map<std::string, std::string> dict;
    std::string input_filename;

};

} // end namespace Utils

#endif /* UTILS_HPP_ */