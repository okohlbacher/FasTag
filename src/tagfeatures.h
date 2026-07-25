#ifndef FASTAG_TAGFEATURES_H
#define FASTAG_TAGFEATURES_H

#include <iosfwd>
#include <string>

namespace FASTag
{
namespace TagFeatures
{

/// Aggregate a FASTag tag TSV into one feature row per spectrum.
bool aggregate(std::istream& input, std::ostream& output, std::string* error = nullptr);

} // namespace TagFeatures
} // namespace FASTag

#endif
