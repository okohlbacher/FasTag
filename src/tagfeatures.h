#ifndef FASTAG_TAGFEATURES_H
#define FASTAG_TAGFEATURES_H

#include <iosfwd>
#include <string>

namespace FasTag
{
namespace TagFeatures
{

/// Aggregate a FasTag tag TSV into one feature row per spectrum.
bool aggregate(std::istream& input, std::ostream& output, std::string* error = nullptr);

} // namespace TagFeatures
} // namespace FasTag

#endif
