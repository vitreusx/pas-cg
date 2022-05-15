#include <cg/input/fields.h>
#include <cg/input/records.h>

namespace cg::records {
std::optional<remark> remark::try_parse(const std::string &line) {
  remark r;

  if (!fields::record_name("REMARK").read(line))
    return std::nullopt;
  r.number = fields::integer(8, 10).read(line);
  r.text = fields::string(12, 80).add(line).read();

  return r;
}

std::string remark::write() const {
  std::string line(80, ' ');

  fields::record_name("REMARK").write(line);
  fields::integer(8, 10).write(line, number);
  fields::string(12, 80).write({&line}, text);

  return line;
}

std::vector<remark> remark::create(int number, const std::string &text) {
  auto max_width = 80 - 12 + 1;
  std::vector<remark> res;
  for (int off = 0; off < (int)text.size(); off += max_width) {
    remark &rem = res.emplace_back();
    rem.number = number;
    rem.text = text.substr(off, max_width);
  }

  return res;
}

std::optional<atom> atom::try_parse(const std::string &line) {
  atom a;

  if (!fields::record_name("ATOM  ").read(line))
    return std::nullopt;
  a.serial = fields::integer(7, 11).read(line);
  a.atom_name = fields::atom(13, 16).read(line);
  a.residue_name = fields::residue_name(18, 20).read(line);
  a.chain_id = fields::achar(22).read(line);
  a.res_seq_num = fields::integer(23, 26).read(line);
  a.pos.x() = fields::real_field(31, 38, 8, 3).read(line);
  a.pos.y() = fields::real_field(39, 46, 8, 3).read(line);
  a.pos.z() = fields::real_field(47, 54, 8, 3).read(line);

  return a;
}

std::string atom::write() const {
  std::string line(80, ' ');

  fields::record_name("ATOM  ").write(line);
  fields::integer(7, 11).write(line, serial);
  fields::atom(13, 16).write(line, atom_name);
  fields::residue_name(18, 20).write(line, residue_name);
  fields::achar(22).write(line, chain_id);
  fields::integer(23, 26).write(line, res_seq_num);
  fields::real_field(31, 38, 8, 3).write(line, pos.x());
  fields::real_field(39, 46, 8, 3).write(line, pos.y());
  fields::real_field(47, 54, 8, 3).write(line, pos.z());

  return line;
}

std::optional<ssbond> ssbond::try_parse(const std::string &line) {
  ssbond s;

  if (!fields::record_name("SSBOND").read(line))
    return std::nullopt;
  s.serial = fields::integer(8, 10).read(line);
  if (!fields::literal(12, 14, "CYS").read(line))
    return std::nullopt;
  s.res[0].chain_id = fields::achar(16).read(line);
  s.res[0].res_seq_num = fields::integer(18, 21).read(line);
  s.res[1].chain_id = fields::achar(30).read(line);
  s.res[1].res_seq_num = fields::integer(32, 35).read(line);
  s.length = fields::real_field(74, 78, 5, 2).read(line);

  return s;
}

std::string ssbond::write() const {
  std::string line(80, ' ');

  fields::record_name("SSBOND").write(line);
  fields::integer(8, 10).write(line, serial);
  fields::literal(12, 14, "CYS").write(line);
  fields::achar(16).write(line, res[0].chain_id);
  fields::integer(18, 21).write(line, res[0].res_seq_num);
  fields::achar(30).write(line, res[1].chain_id);
  fields::integer(32, 35).write(line, res[1].res_seq_num);
  fields::real_field(74, 78, 5, 2).write(line, length);

  return line;
}

std::optional<link> link::try_parse(const std::string &line) {
  link l;

  if (!fields::record_name("LINK  ").read(line))
    return std::nullopt;
  l.res[0].atom_name = fields::atom(13, 16).read(line);
  l.res[0].res_name = fields::residue_name(18, 20).read(line);
  l.res[0].chain_id = fields::achar(22).read(line);
  l.res[0].res_seq_num = fields::integer(23, 26).read(line);
  l.res[1].atom_name = fields::atom(43, 46).read(line);
  l.res[1].res_name = fields::residue_name(48, 50).read(line);
  l.res[1].chain_id = fields::achar(52).read(line);
  l.res[1].res_seq_num = fields::integer(53, 56).read(line);
  l.length = fields::real_field(74, 78, 5, 2).read(line);

  return l;
}

std::string link::write() const {
  std::string line(80, ' ');

  fields::record_name("LINK  ").write(line);
  fields::atom(13, 16).write(line, res[0].atom_name);
  fields::residue_name(18, 20).write(line, res[0].res_name);
  fields::achar(22).write(line, res[0].chain_id);
  fields::integer(23, 26).write(line, res[0].res_seq_num);
  fields::atom(43, 46).write(line, res[1].atom_name);
  fields::residue_name(48, 50).write(line, res[1].res_name);
  fields::achar(52).write(line, res[1].chain_id);
  fields::integer(53, 56).write(line, res[1].res_seq_num);
  fields::real_field(74, 78, 5, 2).write(line, length);

  return line;
}

std::optional<cryst1> cryst1::try_parse(const std::string &line) {
  cryst1 c1;

  if (!fields::record_name("CRYST1").read(line))
    return std::nullopt;
  c1.cell.x() = fields::real_field(7, 15, 9, 3).read(line);
  c1.cell.y() = fields::real_field(16, 24, 9, 3).read(line);
  c1.cell.z() = fields::real_field(25, 33, 9, 3).read(line);

  return c1;
}

std::string cryst1::write() const {
  std::string line(80, ' ');

  fields::record_name("CRYST1").write(line);
  fields::real_field(7, 15, 9, 3).write(line, cell.x());
  fields::real_field(16, 24, 9, 3).write(line, cell.y());
  fields::real_field(25, 33, 9, 3).write(line, cell.z());

  return line;
}

std::optional<ter> ter::try_parse(const std::string &line) {
  ter t;

  if (!fields::record_name("TER   ").read(line))
    return std::nullopt;
  t.serial = fields::integer(7, 11).read(line);
  t.res_name = fields::residue_name(18, 20).read(line);
  t.chain_id = fields::achar(22).read(line);
  t.res_seq_num = fields::integer(23, 26).read(line);

  return t;
}

std::string ter::write() const {
  std::string line(80, ' ');

  fields::record_name("TER   ").write(line);
  fields::integer(7, 11).write(line, serial);
  fields::residue_name(18, 20).write(line, res_name);
  fields::achar(22).write(line, chain_id);
  fields::integer(23, 26).write(line, res_seq_num);

  return line;
}

std::optional<model> model::try_parse(const std::string &line) {
  model m;

  if (!fields::record_name("MODEL ").read(line))
    return std::nullopt;
  m.serial = fields::integer(11, 14).read(line);

  return m;
}

std::string model::write() const {
  std::string line(80, ' ');

  fields::record_name("MODEL ").write(line);
  fields::integer(11, 14).write(line, serial);

  return line;
}

std::optional<endmdl> endmdl::try_parse(const std::string &line) {
  endmdl em;

  if (!fields::record_name("ENDMDL").read(line))
    return std::nullopt;

  return em;
}

std::string endmdl::write() const {
  std::string line(80, ' ');

  fields::record_name("ENDMDL").write(line);

  return line;
}

std::optional<end> end::try_parse(const std::string &line) {
  end e;

  if (!fields::record_name("END   ").read(line))
    return std::nullopt;

  return e;
}

std::string end::write() const {
  std::string line(80, ' ');

  fields::record_name("END   ").write(line);

  return line;
}

std::optional<record> record::try_parse(const std::string &line) {
  if (auto remark_r = remark::try_parse(line); remark_r) {
    return record(remark_r.value());
  } else if (auto atom_r = atom::try_parse(line); atom_r) {
    return record(atom_r.value());
  } else if (auto ssbond_r = ssbond::try_parse(line); ssbond_r) {
    return record(ssbond_r.value());
  } else if (auto link_r = link::try_parse(line); link_r) {
    return record(link_r.value());
  } else if (auto cryst1_r = cryst1::try_parse(line); cryst1_r) {
    return record(cryst1_r.value());
  } else if (auto ter_r = ter::try_parse(line); ter_r) {
    return record(ter_r.value());
  } else if (auto model_r = model::try_parse(line); model_r) {
    return record(model_r.value());
  } else if (auto endmdl_r = endmdl::try_parse(line); endmdl_r) {
    return record(endmdl_r.value());
  } else if (auto end_r = end::try_parse(line); end_r) {
    return record(end_r.value());
  } else {
    return std::nullopt;
  }
}

std::string record::write() const {
  return std::visit(
      [](auto &&value) -> auto { return value.write(); }, rec);
}

record::record(record::record_variant_t rec) : rec{std::move(rec)} {}
} // namespace cg::records