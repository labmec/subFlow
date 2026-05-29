#pragma once

#include "TSFSavable.h"
#include "TPZStream.h"

class TSFFilterCakeMemory : public TSFSavable {
public:
  /// @brief Default constructor
  TSFFilterCakeMemory();

  /// @brief Copy constructor
  TSFFilterCakeMemory(const TSFFilterCakeMemory &other);

  /// @brief  Assignment operator
  /// @param other
  /// @return
  const TSFFilterCakeMemory &operator=(const TSFFilterCakeMemory &other);

  /// @brief  Class identifier
  /// @return
  virtual int ClassId() const;

  /// @brief Destructor
  virtual ~TSFFilterCakeMemory();

  /// @brief Class name
  const std::string Name() const;

  /// @brief Write class attributes
  virtual void Write(TPZStream &buf, int withclassid) const;

  /// @brief Read class attributes
  virtual void Read(TPZStream &buf, void *context);

  /// @brief Print class attributes
  virtual void Print(std::ostream &out = std::cout) const;

  /// @brief Overload of the stream insertion operator
  friend std::ostream &operator<<(std::ostream &out, const TSFFilterCakeMemory &memory) {
    memory.Print(out);
    return out;
  }

  /// @brief  Set the accumulated volume
  /// @param accumulatedVolume
  void SetAccumulatedVolume(REAL accumulatedVolume) { fAccumulatedVolume = accumulatedVolume; }

  /// @brief  Get the accumulated volume
  /// @return
  REAL GetAccumulatedVolume() const { return fAccumulatedVolume; }

protected:
  REAL fAccumulatedVolume;
};