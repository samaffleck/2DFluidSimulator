#pragma once

#include "2DFluidSimulator/Mesh.h"
#include "2DFluidSimulator/IEquation.h"
#include <memory>
#include <vector>

class Domain
{
public:
	Domain() = default;
	~Domain() = default;

	// API
	void addEquation(std::unique_ptr<IEquation> equation);
	void run();
	Mesh& getMesh() { return m_mesh; }
	
private:
	void initialise();

private:
	Mesh m_mesh;
	std::vector<std::unique_ptr<IEquation>> m_equations{};

};
