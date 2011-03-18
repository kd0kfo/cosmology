#include "Help.h"

Help Help::getDefault() {
  Help returnMe;
    returnMe["about"] = "Usage: about x\nDescribes x.";
    returnMe["list"] = "Usage: list\nLists all variables defined.";
    returnMe["define"] = "Usage: defines x y\n Defines x to be the expression y.";
    returnMe["det"] = "Usage: det A\nGives the determinate of A.";
    returnMe["load"] = "Usage: load x <filename>\nLoads the saved Variable as x.";
    returnMe["print"] = "Usage: print x\nPrints x.";
    returnMe["save"] = "Usage: save <filename>\nSaves the previous commands as <filename>.";
    returnMe["remove"] = "Usage: remove x\nRemoves x.";
    returnMe["edit"] = "Usage: edit A m n a,b,...\nEdits the A matrix of size mxn as elements a,b,c,..";
    returnMe["help"] = "Usage: help <topic>\nOffers help on the topic";
    returnMe["last"] = "Usage: last\nPrints the last,  expression used. The last expression is represented by the variable $last$";
    returnMe["dir"] = "Usage: dir [directory]\nlists the contents of the directory.";
    returnMe["quit"] = "Usage: quit\nExits linal.";
    returnMe["exit"] = returnMe["quit"];
    returnMe["physics"] = "Usage: physics <equation> parameter\nGives the output of <equations with the give parameters.";
    returnMe["import"] = "Usage: import\nimports constants.";
    returnMe["set"] = "Usage: set parameter <value>\nSets the parameter to the value <value>.";
    returnMe["run"] = "Usage: run scriptname\n Runs the script of commands which must end with a semicolon.";
    returnMe["functions"] = "List of Functions: ";
    returnMe["history"] = "Lists the past commands.";
    return returnMe;
}

const std::string& Help::getHelp(const arg_t& keyWord) {
    if (this->find(keyWord) != this->end()) {
      return this->at(keyWord);
    } else {
        return "No help for " + keyWord + ".";
    }
}


