#include "processdata.h"

using namespace INMOST;

ProcessData::ProcessData(Mesh *m_)
{
    m = m_;
	rank = m->GetProcessorRank();
}

ProcessData::~ProcessData()
{

}

Mesh *ProcessData::getMesh()
{
    return m;
}
