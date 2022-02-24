#include "fields.h"
#include "utils/text.h"
#include <regex>

namespace cg::fields {
achar::achar(size_t i) : i{i - 1} {}

char achar::read(const std::string &line) const { return line[i]; }

void achar::write(std::string &line, char c) const { line[i] = c; };

atom::atom(size_t i, size_t j) : beg{i - 1}, len{j - i + 1} {}

std::string atom::read(const std::string &line) const {
  return strip(line.substr(beg, len));
}

void atom::write(std::string &line, std::string const &atom) const {
  auto cut_atom = format("%*s", len, atom.c_str());
  format_(line.data() + beg, "%*s", len, atom.c_str());
}

character::character(size_t i) : i{i - 1} {};

char character::read(const std::string &line) const { return line[i]; }

void character::write(std::string &line, char c) const { line[i] = c; }

integer::integer(size_t i, size_t j) : beg{i - 1}, len{j - i + 1} {};

int integer::read(const std::string &line) const {
  return std::stoi(line.substr(beg, len));
}

void integer::write(std::string &line, int n) const {
  auto cut_n = format("%*d", len, n);
  format_(line.data() + beg, "%*d", len, n);
}

lstring::lstring(size_t i, size_t j) : beg{i - 1}, len{j - i + 1} {}

std::string lstring::read(const std::string &line) const {
  return line.substr(beg, len);
}

void lstring::write(std::string &line, std::string const &ls) const {
  format_(line.data() + beg, "%*s", len, ls.c_str());
}

real_field::real_field(size_t i, size_t j, int n, int m)
    : beg{i - 1}, len{j - i + 1}, n{n}, m{m} {}

double real_field::read(const std::string &line) const {
  return std::stod(line.substr(beg, len));
}

void real_field::write(std::string &line, double x) const {
  auto x_str = format("%*.*f", n, m, x);
  format_(line.data() + beg, "%*s", len, x_str.c_str());
}

record_name::record_name(const std::string &name) : name{name} {}

bool record_name::read(const std::string &line) const {
  auto header = lstring(1, 6).read(line);
  return header == name;
}

void record_name::write(std::string &line) const {
  lstring(1, 6).write(line, name);
}

residue_name::residue_name(size_t i, size_t j) : beg{i - 1}, len{j - i + 1} {}

std::string residue_name::read(const std::string &line) const {
  return strip(line.substr(beg, len));
}

void residue_name::write(std::string &line, std::string const &res_name) const {
  format_(line.data() + beg, "%*s", len, res_name.c_str());
}

string::string(size_t i, size_t j) : beg{i - 1}, len{j - i + 1} {}

string &string::add(const std::string &line) {
  cur_text += line;
  return *this;
}

std::string string::read() const {
  return strip(std::regex_replace(cur_text, std::regex("\\s+"), " "));
};

void string::clear() { cur_text = ""; }

void string::write(const std::vector<std::string *> &lines,
                   const std::string &text) const {

  for (size_t line_idx = 0, text_idx = 0; text_idx < text.size();
       ++line_idx, text_idx += len) {
    auto &line = *lines[line_idx];
    format_(line.data() + beg, "%*s", len, &text[text_idx]);
  }
}

literal::literal(size_t i, size_t j, std::string const &text)
    : beg{i - 1}, len{j - i + 1}, text{text} {};

bool literal::read(const std::string &line) const {
  auto line_text = line.substr(beg, len);
  return line_text == text;
}

void literal::write(std::string &line) const {
  format_(line.data() + beg, "%*s", len, text.begin());
}
} // namespace cg::fields