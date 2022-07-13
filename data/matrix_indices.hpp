//
// Created by Bohlender,Ryan J on 7/5/22.
//

#ifndef PERMUTE_ASSOCIATE_MATRIX_INDICES_HPP
#define PERMUTE_ASSOCIATE_MATRIX_INDICES_HPP

enum class Indices : int {
  chrom=0,
  start,
  end,
  ref,
  alt,
  type,
  gene,
  transcript,
  region,
  function,
  annotation,
  first=11
};

#endif // PERMUTE_ASSOCIATE_MATRIX_INDICES_HPP
