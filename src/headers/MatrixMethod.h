//
// Created by mafu on 2024/4/25.
//
#pragma once
// #define EIGEN_NO_DEBUG
// #define EIGEN_USE_MKL_ALL
#include <Eigen/Core>

namespace Coal {
    template<typename ArgsType, typename rowIndexType, typename colIndexType>
    class RefMatrix {
        const ArgsType &args;
        const rowIndexType &row;
        const colIndexType &col;

    public:
        using MatrixType = Eigen::Matrix<
                typename ArgsType::Scalar, rowIndexType::SizeAtCompileTime,
                colIndexType::SizeAtCompileTime,
                ArgsType::Flags & Eigen::RowMajorBit ? Eigen::RowMajor : Eigen::ColMajor,
                rowIndexType::MaxSizeAtCompileTime, colIndexType::MaxSizeAtCompileTime>;
        RefMatrix(const ArgsType &args, const rowIndexType &row, const colIndexType &col)
            : args(args), row(row), col(col) {}

        const typename ArgsType::Scalar &operator()(Eigen::Index i,
                                                    Eigen::Index j) const {
            return args(row[i], col[j]);
        }
    };
    template<typename ArgsType, typename rowIndexType, typename colIndexType>
    Eigen::CwiseNullaryOp<
            RefMatrix<ArgsType, rowIndexType, colIndexType>,
            typename RefMatrix<ArgsType, rowIndexType, colIndexType>::MatrixType>
    selectByIndex(const Eigen::MatrixBase<ArgsType> &arg, const rowIndexType &row,
                  const colIndexType &col) {
        using Func       = RefMatrix<ArgsType, rowIndexType, colIndexType>;
        using MatrixType = typename Func::MatrixType;
        return MatrixType::NullaryExpr(row.size(), col.size(),
                                       Func(arg.derived(), row, col));
    }
    template<typename Derived>
    Eigen::ArrayXi findRowsGreaterThanV2(const Eigen::MatrixBase<Derived> &matrix,
                                         int colIndex,
                                         const typename Derived::Scalar &value) {
        Eigen::ArrayXi rowIndexes(matrix.rows());
        int count = 0;
        if (colIndex < 0 || colIndex >= matrix.cols()) {
            throw std::out_of_range("Column index is out of range.");
        }
        for (int i = 0; i < matrix.rows(); ++i) {
            if (matrix(i, colIndex) > value) {
                rowIndexes(count++) = i;
            }
        }
        rowIndexes.conservativeResize(count);
        return rowIndexes;
    }
    template<typename Derived, typename OtherDerived>
    void conditionSelect(const Eigen::MatrixBase<Derived> &tempMatrixs,
                         const typename Derived::Scalar &prob_cut,
                         Eigen::MatrixBase<OtherDerived> &matrixSwap) {
        const Eigen::ArrayXi indices = findRowsGreaterThanV2(tempMatrixs, 10, prob_cut);
        const Eigen::ArrayXi indices_col = Eigen::ArrayXi::LinSpaced(11, 0, 10);
        matrixSwap.derived().resize(indices.size(), 11);
        matrixSwap.derived() = selectByIndex(tempMatrixs, indices, indices_col);
    }

    template<typename Derived, typename OtherDerived>
    void select(const Eigen::MatrixBase<Derived> &matrix, int colIndex,
                typename Derived::Scalar value,
                Eigen::MatrixBase<OtherDerived> &matrixSwap) {

        Eigen::ArrayXi rowIndices(matrix.rows());
        int count = 0;
        for (int i = 0; i < matrix.rows(); ++i) {
            if (matrix(i, colIndex) > value) {
                rowIndices(count++) = i;
            }
        }
        rowIndices.conservativeResize(count);
        const Eigen::ArrayXi colIndices =
                Eigen::ArrayXi::LinSpaced(matrix.cols(), 0, matrix.cols() - 1);

        matrixSwap.derived().resize(rowIndices.size(), matrix.cols());
        matrixSwap.derived().noalias() = selectByIndex(matrix, rowIndices, colIndices);
    }

    template<typename MatrixTypeA, typename MatrixTypeB>
    void appendMatrixToBottom(const Eigen::MatrixBase<MatrixTypeA> &src,
                              Eigen::MatrixBase<MatrixTypeB> &dest, long &numRowsInDest) {
        int newRows = numRowsInDest + src.rows();
        if (newRows > dest.rows()) {
            constexpr int addSize = 100000;
            dest.derived().conservativeResize(newRows + addSize, Eigen::NoChange);
        }
        dest.block(numRowsInDest, 0, src.rows(), src.cols()) = src;
        numRowsInDest += src.rows();
    }

}// namespace Coal
