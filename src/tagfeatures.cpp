// tagfeatures converts FasTag's one-row-per-tag TSV output into a deterministic
// one-row-per-spectrum feature table for joining tag evidence to downstream PSM
// rescoring inputs; it is intentionally standalone and uses only the C++ standard
// library.

#include "tagfeatures.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{

std::vector<std::string> splitTabs(const std::string& line)
{
  std::vector<std::string> fields;
  std::size_t begin = 0;
  for (;;)
  {
    const std::size_t tab = line.find('\t', begin);
    if (tab == std::string::npos)
    {
      fields.push_back(line.substr(begin));
      return fields;
    }
    fields.push_back(line.substr(begin, tab - begin));
    begin = tab + 1;
  }
}

bool parseDouble(const std::string& text, double& value)
{
  if (text.empty()) return false;
  char* end = nullptr;
  errno = 0;
  value = std::strtod(text.c_str(), &end);
  return end != text.c_str() && *end == '\0' && std::isfinite(value);
}

bool parseSize(const std::string& text, std::size_t& value)
{
  if (text.empty() || text.front() == '-') return false;
  char* end = nullptr;
  errno = 0;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE ||
      parsed > std::numeric_limits<std::size_t>::max())
  {
    return false;
  }
  value = static_cast<std::size_t>(parsed);
  return true;
}

bool parseFlag(const std::string& text, bool& value)
{
  if (text == "0")
  {
    value = false;
    return true;
  }
  if (text == "1")
  {
    value = true;
    return true;
  }
  return false;
}

struct Columns
{
  std::size_t spectrum = 0;
  std::size_t length = 0;
  std::size_t evalue = 0;
  std::size_t gapped = 0;
  std::size_t extended = 0;
  std::size_t min_conf = 0;
  std::size_t mean_conf = 0;
  std::size_t fasta_hit = 0;
  bool has_min_conf = false;
  bool has_mean_conf = false;
  bool has_fasta_hit = false;
};

bool findColumns(const std::vector<std::string>& header, Columns& columns,
                 std::string* error)
{
  std::unordered_map<std::string, std::size_t> index;
  for (std::size_t i = 0; i < header.size(); ++i) index.emplace(header[i], i);

  auto required = [&](const char* name, std::size_t& destination)
  {
    const auto found = index.find(name);
    if (found == index.end())
    {
      if (error) *error = std::string("missing required column '") + name + "'";
      return false;
    }
    destination = found->second;
    return true;
  };

  if (!required("spectrum", columns.spectrum) ||
      !required("length", columns.length) ||
      !required("evalue", columns.evalue) ||
      !required("gapped", columns.gapped) ||
      !required("extended", columns.extended))
  {
    return false;
  }

  auto optional = [&](const char* name, std::size_t& destination)
  {
    const auto found = index.find(name);
    if (found == index.end()) return false;
    destination = found->second;
    return true;
  };

  columns.has_min_conf = optional("min_conf", columns.min_conf);
  columns.has_mean_conf = optional("mean_conf", columns.mean_conf);
  columns.has_fasta_hit = optional("fasta_hit", columns.fasta_hit);
  return true;
}

struct Aggregate
{
  std::size_t n_tags = 0;
  double best_evalue = std::numeric_limits<double>::infinity();
  double second_best_evalue = std::numeric_limits<double>::infinity();
  std::size_t longest_tag = 0;
  double best_min_conf = -std::numeric_limits<double>::infinity();
  long double sum_mean_conf = 0.0L;
  std::size_t n_gapped = 0;
  std::size_t n_extended = 0;
  std::size_t n_fasta_hits = 0;

  void add(double evalue, std::size_t length, bool gapped, bool extended,
           double min_conf, double mean_conf, bool fasta_hit)
  {
    ++n_tags;
    if (evalue < best_evalue)
    {
      second_best_evalue = best_evalue;
      best_evalue = evalue;
    }
    else if (evalue < second_best_evalue)
    {
      second_best_evalue = evalue;
    }
    longest_tag = std::max(longest_tag, length);
    best_min_conf = std::max(best_min_conf, min_conf);
    sum_mean_conf += mean_conf;
    n_gapped += gapped ? 1 : 0;
    n_extended += extended ? 1 : 0;
    n_fasta_hits += fasta_hit ? 1 : 0;
  }
};

} // namespace

namespace FasTag
{
namespace TagFeatures
{

bool aggregate(std::istream& input, std::ostream& output, std::string* error)
{
  std::string line;
  if (!std::getline(input, line))
  {
    if (error) *error = "input is empty";
    return false;
  }
  if (!line.empty() && line.back() == '\r') line.pop_back();
  if (line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);

  Columns columns;
  if (!findColumns(splitTabs(line), columns, error)) return false;

  std::map<std::string, Aggregate> spectra;
  while (std::getline(input, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    const std::vector<std::string> fields = splitTabs(line);

    std::size_t length = 0;
    double evalue = 0.0;
    bool gapped = false;
    bool extended = false;
    double min_conf = 0.0;
    double mean_conf = 0.0;
    bool fasta_hit = false;

    auto has = [&](std::size_t column) { return column < fields.size(); };
    if (!has(columns.spectrum) || fields[columns.spectrum].empty() ||
        !has(columns.length) || !parseSize(fields[columns.length], length) ||
        !has(columns.evalue) || !parseDouble(fields[columns.evalue], evalue) ||
        evalue < 0.0 ||
        !has(columns.gapped) || !parseFlag(fields[columns.gapped], gapped) ||
        !has(columns.extended) || !parseFlag(fields[columns.extended], extended))
    {
      continue;
    }
    if (columns.has_min_conf &&
        (!has(columns.min_conf) || !parseDouble(fields[columns.min_conf], min_conf)))
    {
      continue;
    }
    if (columns.has_mean_conf &&
        (!has(columns.mean_conf) || !parseDouble(fields[columns.mean_conf], mean_conf)))
    {
      continue;
    }
    if (columns.has_fasta_hit)
    {
      if (!has(columns.fasta_hit) || fields[columns.fasta_hit].empty()) continue;
      fasta_hit = fields[columns.fasta_hit] != "-";
    }

    spectra[fields[columns.spectrum]].add(
        evalue, length, gapped, extended, min_conf, mean_conf, fasta_hit);
  }

  output << "spectrum\tn_tags\tbest_evalue\tsecond_best_evalue\tevalue_gap"
            "\tlongest_tag\tbest_min_conf\tmean_mean_conf\tfrac_gapped"
            "\tfrac_extended\tn_fasta_hits\n";
  output << std::setprecision(17);
  for (const auto& item : spectra)
  {
    const Aggregate& value = item.second;
    const double second = value.n_tags == 1
                              ? value.best_evalue
                              : value.second_best_evalue;
    const double gap = std::log10(second + 1e-300) -
                       std::log10(value.best_evalue + 1e-300);
    output << item.first << '\t'
           << value.n_tags << '\t'
           << value.best_evalue << '\t'
           << second << '\t'
           << gap << '\t'
           << value.longest_tag << '\t';
    if (columns.has_min_conf) output << value.best_min_conf;
    else output << "NA";
    output << '\t';
    if (columns.has_mean_conf)
      output << value.sum_mean_conf / static_cast<long double>(value.n_tags);
    else
      output << "NA";
    output << '\t'
           << static_cast<double>(value.n_gapped) / value.n_tags << '\t'
           << static_cast<double>(value.n_extended) / value.n_tags << '\t'
           << value.n_fasta_hits << '\n';
  }

  if (!output)
  {
    if (error) *error = "failed while writing output";
    return false;
  }
  return true;
}

} // namespace TagFeatures
} // namespace FasTag

#ifndef TAGFEATURES_NO_MAIN
int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "usage: tagfeatures <fastag-tags.tsv> <features.tsv>\n";
    return 2;
  }

  std::ifstream input(argv[1]);
  if (!input)
  {
    std::cerr << "tagfeatures: cannot open input '" << argv[1] << "'\n";
    return 1;
  }
  std::ofstream output(argv[2]);
  if (!output)
  {
    std::cerr << "tagfeatures: cannot open output '" << argv[2] << "'\n";
    return 1;
  }

  std::string error;
  if (!FasTag::TagFeatures::aggregate(input, output, &error))
  {
    std::cerr << "tagfeatures: " << error << '\n';
    return 1;
  }
  return 0;
}
#endif
