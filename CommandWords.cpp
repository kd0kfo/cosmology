#include "CommandWords.h"

bool CommandWords::isCommand(const arg_t& aString)
{
    size_t size = this->size();
    for (size_t i = 0; i < size; i++)
        if (this->at(i) == aString)
            return true;
    return false;
}

arg_t CommandWords::showAll() {
    arg_t bean = "";

    std::set<std::string> sorter;
    args_t::const_iterator words = this->begin();
    for(;words != this->end();words++)
      sorter.insert(*words);
    std::set<std::string>::const_iterator word = sorter.begin();
    for(;word != sorter.end();word++)
      bean += " " + *word;

    return bean + " and functions";
}
