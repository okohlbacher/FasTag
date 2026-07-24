#include "tagfeatures.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static int failures = 0;
#define CHECK(c, ...) do { if (!(c)) { ++failures; std::printf("FAIL %d: ", __LINE__); \
  std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static std::vector<std::string> splitTabs(const std::string& line)
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

static bool near(double actual, double expected, double tolerance = 1e-12)
{
  return std::fabs(actual - expected) <= tolerance;
}

static std::vector<std::vector<std::string>> aggregateFile(
    const std::filesystem::path& path, std::string& header)
{
  std::ifstream input(path);
  std::ostringstream output;
  std::string error;
  const bool ok = FasTag::TagFeatures::aggregate(input, output, &error);
  CHECK(ok, "aggregation failed: %s", error.c_str());

  std::istringstream result(output.str());
  std::getline(result, header);
  std::vector<std::vector<std::string>> rows;
  std::string line;
  while (std::getline(result, line)) rows.push_back(splitTabs(line));
  return rows;
}

int main()
{
  const auto unique = std::chrono::high_resolution_clock::now()
                          .time_since_epoch().count();
  const std::filesystem::path full_path =
      std::filesystem::temp_directory_path() /
      ("tagfeatures_test_" + std::to_string(unique) + ".tsv");
  {
    std::ofstream out(full_path);
    out << "spectrum\ttag\tlength\tcharge\tnterm_mass\tcterm_mass\textended"
           "\tgapped\tevalue\tmin_conf\tmean_conf\tfasta_hit\n"
        << "z_spec\tZZZ\t4\t2\t0\t0\t0\t0\t0.5\t0.2\t0.4\t-\n"
        << "a_spec\tAAA\t3\t2\t0\t0\t1\t0\t1e-5\t0.4\t0.6\t-\n"
        << "a_spec\tBBBBB\t5\t2\t0\t0\t0\t1\t1e-2\t0.8\t0.9\tP12345\n"
        << "z_spec\tZZZZZ\t5\t2\t0\t0\t1\t1\t0.7\t0.3\t0.5\tZ9\n"
        << "a_spec\tBAD\t7\t2\t0\t0\t1\t1\tnot-a-number\t1\t1\tBAD\n"
        << "a_spec\tCCCC\t4\t2\t0\t0\t1\t1\t1e-3\t0.6\t0.3\tQ67890\n"
        << "malformed\ttoo\tshort\n";
  }

  std::string header;
  const auto rows = aggregateFile(full_path, header);
  CHECK(header ==
            "spectrum\tn_tags\tbest_evalue\tsecond_best_evalue\tevalue_gap"
            "\tlongest_tag\tbest_min_conf\tmean_mean_conf\tfrac_gapped"
            "\tfrac_extended\tn_fasta_hits",
        "unexpected output header");
  CHECK(rows.size() == 2, "expected 2 spectra, got %zu", rows.size());
  if (rows.size() >= 2 && rows[0].size() == 11)
  {
    const auto& a = rows[0];
    CHECK(a[0] == "a_spec", "spectra are not lexicographically ordered");
    CHECK(a[1] == "3", "n_tags: expected 3, got %s", a[1].c_str());
    CHECK(near(std::stod(a[2]), 1e-5), "best_evalue mismatch");
    CHECK(near(std::stod(a[3]), 1e-3), "second_best_evalue mismatch");
    CHECK(near(std::stod(a[4]), 2.0), "evalue_gap mismatch");
    CHECK(a[5] == "5", "longest_tag: expected 5, got %s", a[5].c_str());
    CHECK(near(std::stod(a[6]), 0.8), "best_min_conf mismatch");
    CHECK(near(std::stod(a[7]), 0.6), "mean_mean_conf mismatch");
    CHECK(near(std::stod(a[8]), 2.0 / 3.0), "frac_gapped mismatch");
    CHECK(near(std::stod(a[9]), 2.0 / 3.0), "frac_extended mismatch");
    CHECK(a[10] == "2", "n_fasta_hits: expected 2, got %s", a[10].c_str());
  }
  else if (!rows.empty())
  {
    CHECK(false, "first feature row has %zu columns", rows[0].size());
  }

  const std::filesystem::path old_path =
      std::filesystem::temp_directory_path() /
      ("tagfeatures_old_test_" + std::to_string(unique) + ".tsv");
  {
    std::ofstream out(old_path);
    out << "spectrum\ttag\tlength\tcharge\tnterm_mass\tcterm_mass\textended"
           "\tgapped\tevalue\n"
        << "old_spec\tOLD\t3\t2\t0\t0\t0\t1\t0.1\n";
  }

  std::string old_header;
  const auto old_rows = aggregateFile(old_path, old_header);
  CHECK(old_rows.size() == 1, "expected one old-format spectrum");
  if (old_rows.size() == 1 && old_rows[0].size() == 11)
  {
    const auto& old = old_rows[0];
    CHECK(old[3] == old[2], "single tag must reuse best as second best");
    CHECK(near(std::stod(old[4]), 0.0), "single-tag evalue_gap must be zero");
    CHECK(old[6] == "NA", "missing min_conf must produce NA");
    CHECK(old[7] == "NA", "missing mean_conf must produce NA");
    CHECK(near(std::stod(old[8]), 1.0), "old-format frac_gapped mismatch");
    CHECK(near(std::stod(old[9]), 0.0), "old-format frac_extended mismatch");
    CHECK(old[10] == "0", "missing fasta_hit must produce zero hits");
  }

  std::error_code ignored;
  std::filesystem::remove(full_path, ignored);
  std::filesystem::remove(old_path, ignored);

  if (failures)
  {
    std::printf("tagfeatures_test: %d checks failed\n", failures);
    return 1;
  }
  std::printf("tagfeatures_test: all checks passed\n");
  return 0;
}
