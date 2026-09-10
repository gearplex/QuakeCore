#include "quake/update_basis.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace quake {

SparseUpdateBasis::SparseUpdateBasis(int rows, int cols, std::vector<int> col_ptr,
                                     std::vector<int> row_ind, std::vector<double> values)
    : rows_(rows), cols_(cols), col_ptr_(std::move(col_ptr)),
      row_ind_(std::move(row_ind)), values_(std::move(values)) {
    if (rows_ < 0 || cols_ < 0 || static_cast<int>(col_ptr_.size()) != cols_ + 1 ||
        row_ind_.size() != values_.size() || col_ptr_.empty() || col_ptr_.front() != 0 ||
        col_ptr_.back() != static_cast<int>(values_.size()))
        throw std::invalid_argument("SparseUpdateBasis dimension mismatch");
    for (int c=0;c<cols_;++c) {
        if (col_ptr_[static_cast<std::size_t>(c)] > col_ptr_[static_cast<std::size_t>(c+1)])
            throw std::invalid_argument("SparseUpdateBasis invalid column pointers");
        int prev=-1;
        for(int p=col_ptr_[static_cast<std::size_t>(c)];p<col_ptr_[static_cast<std::size_t>(c+1)];++p){
            const int r=row_ind_[static_cast<std::size_t>(p)];
            if(r<0||r>=rows_||r<=prev) throw std::invalid_argument("SparseUpdateBasis row indices must be sorted/unique");
            if(!std::isfinite(values_[static_cast<std::size_t>(p)])) throw std::invalid_argument("SparseUpdateBasis nonfinite value");
            prev=r;
        }
    }
}

SparseUpdateBasis SparseUpdateBasis::from_dense(int rows,int cols,
                                                 const std::vector<double>& a,double tol){
    if(rows<0||cols<0||static_cast<int>(a.size())!=rows*cols||tol<0.0)
        throw std::invalid_argument("SparseUpdateBasis dense input mismatch");
    std::vector<std::vector<std::pair<int,double>>> columns(static_cast<std::size_t>(cols));
    for(int c=0;c<cols;++c) for(int r=0;r<rows;++r){
        const double v=a[static_cast<std::size_t>(c*rows+r)];
        if(std::abs(v)>tol) columns[static_cast<std::size_t>(c)].push_back({r,v});
    }
    return from_columns(rows,columns,tol);
}

SparseUpdateBasis SparseUpdateBasis::from_columns(
    int rows,const std::vector<std::vector<std::pair<int,double>>>& columns,double tol){
    if(rows<0||tol<0.0) throw std::invalid_argument("SparseUpdateBasis invalid input");
    std::vector<int> cp(columns.size()+1,0),ri; std::vector<double> vv;
    for(std::size_t c=0;c<columns.size();++c){
        std::map<int,double> merged;
        for(const auto& [r,v]:columns[c]){
            if(r<0||r>=rows||!std::isfinite(v)) throw std::invalid_argument("SparseUpdateBasis column entry");
            merged[r]+=v;
        }
        for(const auto& [r,v]:merged) if(std::abs(v)>tol){ri.push_back(r);vv.push_back(v);}
        cp[c+1]=static_cast<int>(vv.size());
    }
    return SparseUpdateBasis(rows,static_cast<int>(columns.size()),std::move(cp),std::move(ri),std::move(vv));
}

double SparseUpdateBasis::column_dot(int col,const std::vector<double>& x) const{
    if(col<0||col>=cols_||static_cast<int>(x.size())!=rows_) throw std::invalid_argument("SparseUpdateBasis dot dimension");
    double s=0.0;for(int p=col_ptr_[static_cast<std::size_t>(col)];p<col_ptr_[static_cast<std::size_t>(col+1)];++p)
        s+=values_[static_cast<std::size_t>(p)]*x[static_cast<std::size_t>(row_ind_[static_cast<std::size_t>(p)])];
    return s;
}

std::vector<double> SparseUpdateBasis::dense_column(int col) const{
    if(col<0||col>=cols_) throw std::out_of_range("SparseUpdateBasis column");
    std::vector<double> out(static_cast<std::size_t>(rows_),0.0);
    axpy_column(col,1.0,out);return out;
}

std::vector<double> SparseUpdateBasis::dense_column_major() const{
    std::vector<double> out(static_cast<std::size_t>(rows_*cols_),0.0);
    for(int c=0;c<cols_;++c) for(int p=col_ptr_[static_cast<std::size_t>(c)];p<col_ptr_[static_cast<std::size_t>(c+1)];++p)
        out[static_cast<std::size_t>(c*rows_+row_ind_[static_cast<std::size_t>(p)])]=values_[static_cast<std::size_t>(p)];
    return out;
}

void SparseUpdateBasis::axpy_column(int col,double scale,std::vector<double>& y) const{
    if(col<0||col>=cols_||static_cast<int>(y.size())!=rows_) throw std::invalid_argument("SparseUpdateBasis axpy dimension");
    for(int p=col_ptr_[static_cast<std::size_t>(col)];p<col_ptr_[static_cast<std::size_t>(col+1)];++p)
        y[static_cast<std::size_t>(row_ind_[static_cast<std::size_t>(p)])]+=scale*values_[static_cast<std::size_t>(p)];
}

double SparseUpdateBasis::value(int row,int col) const{
    if(row<0||row>=rows_||col<0||col>=cols_) throw std::out_of_range("SparseUpdateBasis value");
    const int b=col_ptr_[static_cast<std::size_t>(col)],e=col_ptr_[static_cast<std::size_t>(col+1)];
    auto it=std::lower_bound(row_ind_.begin()+b,row_ind_.begin()+e,row);
    if(it==row_ind_.begin()+e||*it!=row) return 0.0;
    return values_[static_cast<std::size_t>(std::distance(row_ind_.begin(),it))];
}

} // namespace quake
