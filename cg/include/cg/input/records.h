#pragma once
#include <cg/types/vec3.h>
#include <iostream>
#include <optional>
#include <variant>

namespace cg::records {
class remark {
public:
  int number;
  std::string text;

public:
  static std::optional<remark> try_parse(std::string const &line);
  std::string write() const;

  static std::vector<remark> create(int number, std::string const &text);
};

class atom {
public:
  std::string serial;
  std::string atom_name;
  std::string residue_name;
  std::string chain_id;
  int res_seq_num;
  vec3<double> pos;

public:
  static std::optional<atom> try_parse(std::string const &line);
  std::string write() const;
};

class ssbond {
public:
  std::string serial;
  struct {
    std::string chain_id;
    int res_seq_num;
  } res[2];
  double length;

public:
  static std::optional<ssbond> try_parse(std::string const &line);
  std::string write() const;
};

class link {
public:
  struct {
    std::string atom_name;
    std::string res_name;
    std::string chain_id;
    int res_seq_num;
  } res[2];
  double length;

public:
  static std::optional<link> try_parse(std::string const &line);
  std::string write() const;
};

class cryst1 {
public:
  vec3<double> cell;
  double alpha = 90.0, beta = 90.0, gamma = 90.0;
  std::string sgroup = "P 1";
  int z = 1;

public:
  static std::optional<cryst1> try_parse(std::string const &line);
  std::string write() const;
};

class ter {
public:
  std::string serial;
  std::string res_name;
  std::string chain_id;
  int res_seq_num;

public:
  static std::optional<ter> try_parse(std::string const &line);
  std::string write() const;
};

class model {
public:
  int serial;

public:
  static std::optional<model> try_parse(std::string const &line);
  std::string write() const;
};

class endmdl {
public:
  static std::optional<endmdl> try_parse(std::string const &line);
  std::string write() const;
};

class end {
public:
  static std::optional<end> try_parse(std::string const &line);
  std::string write() const;
};

class record {
public:
  using record_variant_t =
      std::variant<remark, atom, ssbond, link, cryst1, ter, model, endmdl, end>;
  record_variant_t rec;

public:
  static std::optional<record> try_parse(std::string const &line);
  std::string write() const;

  template <typename Record>
  Record *cast() {
    return std::get_if<Record>(&rec);
  }

  template <typename Record>
  Record const *cast() const {
    return std::get_if<Record>(&rec);
  }

private:
  explicit record(record_variant_t rec);
};
} // namespace cg::records
