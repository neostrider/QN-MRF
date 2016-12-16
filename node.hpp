#ifndef NODE_HPP
#define NODE_HPP

class node
{
    int nLabel_;
    std::vector<double> primalFrac_;

public:

    node(int nLabel):nLabel_(nLabel) {
    }

    int getLabelSiz() {
     return nLabel_;
    }

    std::vector<double> getPrimalFrac() {
     return primalFrac_;
    }
};

#endif // NODE_HPP
