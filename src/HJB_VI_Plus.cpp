#include "HJB_Stability.h"

const double PI = 3.141592653589793238463;

// The lower and upper bound of the pendulum variables
double LLow = 0.35;             double LUpp = 1.05;
double LdotLow = -3.0;          double LdotUpp = 3.0;
double ThetaLow = -PI/2.0;      double ThetaUpp = -1.0 * ThetaLow;
double ThetadotLow = -3.0;      double ThetadotUpp = -1.0 * ThetadotLow;

// const int L_Grids = 71;
const int L_Grids = 61;
const int Ldot_Grids = 61;
// const int Theta_Grids = 181;
// const int Thetadot_Grids = 121;
const int Theta_Grids = 61;
const int Thetadot_Grids = 61;

double Length_Low = 0.35;
double Length_Upp = 1.05;

double g;

const int F_Grids = 51;       // This is the discretization of the control points within a range
double delta_t = 0.05;

float Ctv = 10.0;     // coefficient of theta violation
float Clv = 10.0;     // coefficient of length violation

double Clow = 0.5;
double Cupp = 2.0;


// The corresponding dimension size of the StateMatrix
double L_length = LUpp - LLow;
double Ldot_length = LdotUpp - LdotLow;
double Theta_length = ThetaUpp - ThetaLow;
double Thetadot_length = ThetadotUpp - ThetadotLow;

double L_unit = L_length/(1.0* L_Grids - 1.0);
double Ldot_unit = Ldot_length/(1.0* Ldot_Grids - 1.0);
double Theta_unit = Theta_length/(1.0* Theta_Grids - 1.0);
double Thetadot_unit = Thetadot_length/(1.0* Thetadot_Grids - 1.0);

std::vector<double> L_vector(L_Grids), Ldot_vector(Ldot_Grids), Theta_vector(Theta_Grids), Thetadot_vector(Thetadot_Grids);

void StateVectorSubs(std::vector<double>& L_vector,std::vector<double>& Ldot_vector,std::vector<double>& Theta_vector, std::vector<double>& Thetadot_vector)
{
  // L_vector filling
  for (int i = 0; i < L_Grids; i++)
  {
    L_vector[i] = LLow + (1.0 * i) * L_unit;
  }
  // Ldot_vector filling
  for (int i = 0; i < Ldot_Grids; i++)
  {
    Ldot_vector[i] = LdotLow + (1.0 * i) * Ldot_unit;
  }
  // Theta_vector filling
  for (int i = 0; i < Theta_Grids; i++)
  {
    Theta_vector[i] = ThetaLow + (1.0 * i) * Theta_unit;
  }
  // Thetadot_vector filling
  for (int i = 0; i < Thetadot_Grids; i++)
  {
    Thetadot_vector[i] = ThetadotLow + (1.0 * i) * Thetadot_unit;
  }
}
SystemIndex SystemState2Index(const SystemState & x)
{
  // This function is used to transform the System State from double into the corresponding index
  // Here the x must have already been checked with the given pendulum length bound.

  double L_FloatIndex =         (x.L - LLow)/L_unit * 1.0;
  double Ldot_FloatIndex =      (x.Ldot - LdotLow)/Ldot_unit * 1.0;
  double Theta_FloatIndex =     (x.Theta - ThetaLow)/Theta_unit * 1.0;
  double Thetadot_FloatIndex =  (x.Thetadot - ThetadotLow)/Thetadot_unit * 1.0;

  int L_Index =         std::round(L_FloatIndex);
  int Ldot_Index =      std::round(Ldot_FloatIndex);
  int Theta_Index =     std::round(Theta_FloatIndex);
  int Thetadot_Index =  std::round(Thetadot_FloatIndex);

  // The L Index reshift
  if(L_Index<0)
  {
    L_Index = 0;
  }
  if(L_Index>=L_Grids)
  {
    L_Index = L_Grids-1;
  }

  // The Ldot Index reshift
  if(Ldot_Index<0)
  {
    Ldot_Index = 0;
  }
  if(Ldot_Index>=Ldot_Grids)
  {
    Ldot_Index = Ldot_Grids-1;
  }

  // The Theta Index reshift
  if(Theta_Index<0)
  {
    Theta_Index = 0;
  }
  if(Theta_Index>=Theta_Grids)
  {
    Theta_Index = Theta_Grids-1;
  }

  // The Thetadot Index reshift
  if(Thetadot_Index<0)
  {
    Thetadot_Index = 0;
  }
  if(Thetadot_Index>=Thetadot_Grids)
  {
    Thetadot_Index = Thetadot_Grids-1;
  }

  SystemIndex SystemIndex_i(L_Index, Ldot_Index, Theta_Index, Thetadot_Index);
  return SystemIndex_i;
}
int LdotValue2Index(const double& Ldot)
{
  double Ldot_FloatIndex =      (Ldot - LdotLow)/Ldot_unit * 1.0;
  int Ldot_Index =              std::round(Ldot_FloatIndex);
  // The Ldot Index reshift
  if(Ldot_Index<0)
  {
    Ldot_Index = 0;
  }
  if(Ldot_Index>=Ldot_Grids)
  {
    Ldot_Index = Ldot_Grids-1;
  }
  return Ldot_Index;
}
SystemIndex SystemState2StateIndex(const SystemState & x)
{
  // This function is used to transform the System State from double into the corresponding index
  // Here the x must have already been checked with the given pendulum length bound.

  double L_FloatIndex =         (x.L - LLow)/L_unit * 1.0;
  double Ldot_FloatIndex =      (x.Ldot - LdotLow)/Ldot_unit * 1.0;
  double Theta_FloatIndex =     (x.Theta - ThetaLow)/Theta_unit * 1.0;
  double Thetadot_FloatIndex =  (x.Thetadot - ThetadotLow)/Thetadot_unit * 1.0;

  int L_Index =         std::round(L_FloatIndex);
  int Ldot_Index =      std::round(Ldot_FloatIndex);
  int Theta_Index =     std::round(Theta_FloatIndex);
  int Thetadot_Index =  std::round(Thetadot_FloatIndex);

  // The L Index reshift
  if(L_Index<0)
  {
    L_Index = 0;
  }
  if(L_Index>=L_Grids)
  {
    L_Index = L_Grids-1;
  }

  // The Ldot Index reshift
  if(Ldot_Index<0)
  {
    Ldot_Index = 0;
  }
  if(Ldot_Index>=Ldot_Grids)
  {
    Ldot_Index = Ldot_Grids-1;
  }

  // The Theta Index reshift
  if(Theta_Index<0)
  {
    Theta_Index = 0;
  }
  if(Theta_Index>=Theta_Grids)
  {
    Theta_Index = Theta_Grids-1;
  }

  // The Thetadot Index reshift
  if(Thetadot_Index<0)
  {
    Thetadot_Index = 0;
  }
  if(Thetadot_Index>=Thetadot_Grids)
  {
    Thetadot_Index = Thetadot_Grids-1;
  }

  SystemIndex SystemIndex_i(L_Index, Ldot_Index, Theta_Index, Thetadot_Index);
  return SystemIndex_i;
}
std::tuple<int, int, int, int> SystemState2Tuple(const SystemState & x, const std::vector<double> & L_vector,const std::vector<double> & Ldot_vector,const std::vector<double> & Theta_vector, const std::vector<double> & Thetadot_vector)
{
  SystemIndex x_index = SystemState2StateIndex(x);
  std::tuple<int, int, int, int> Tuple_k(x_index.L_index, x_index.Ldot_index, x_index.Theta_index, x_index.Thetadot_index);
  return Tuple_k;
}
std::tuple<int, int, int, int> SystemIndex2Tuple(const SystemIndex& x_index)
{
  std::tuple<int, int, int, int> Tuple_k(x_index.L_index, x_index.Ldot_index, x_index.Theta_index, x_index.Thetadot_index);
  return Tuple_k;
}
void CostUpdate(DBNode& Node_i)
{
  // KE obj
  double L_i, Ldot_i, Theta_i, Thetadot_i;
  L_i = Node_i.NodeState.L;
  Ldot_i = Node_i.NodeState.Ldot;
  Theta_i = Node_i.NodeState.Theta;
  Thetadot_i = Node_i.NodeState.Thetadot;

  if((Theta_i>=0)&&(Thetadot_i>0))
  {
    Node_i.cost = 0.0;
    return;
  }

  float KE_ijkl, lv_ijkl, tv_ijkl;

  KE_ijkl = Ldot_i * Ldot_i  + L_i * L_i * Thetadot_i * Thetadot_i;

  // length violation obj
  double Node_k_L = Node_i.NodeState.L - Node_i.NodeState.Ldot * delta_t;

  switch (Node_i.LengthFeasibleFlag)
  {
    case 0:
    // Length Infeasible Cost
    if(Node_k_L<LLow)
    {
      lv_ijkl = Clv * (LLow - Node_k_L);
    }
    else
    {
      lv_ijkl = Clv * (Node_k_L - LUpp);
    }
    break;
    case 1:
    lv_ijkl = 0.0;
    break;
  }
  // theta violation obj
  switch (Node_i.ThetaViableFlag)
  {
    case 0:
    tv_ijkl = -1.0 * Ctv * Theta_i;
    break;
    case 1:
    tv_ijkl = 0.0;
    break;
  }
  // Total objective function
  Node_i.cost = KE_ijkl + lv_ijkl + tv_ijkl;
}
void DBNodeTLcUpdate(DBNode & Node_i)
{
  // This function can only be called after Node_i has been initialized
  // The following properties for Node_i will be updated.
  //  1. ThetaViableFlag
  //  2. LengthFeasibleFlag
  //  3. cost

  // 1. ThetaViableFlag
  if(Node_i.NodeState.Theta<0){  Node_i.ThetaViableFlag = 0;  }
  // 2. LengthFeasibleFlag
  double Node_k_L = Node_i.NodeState.L - Node_i.NodeState.Ldot * delta_t;
  if((Node_k_L<LLow)||(Node_k_L>LUpp))
  {
    Node_i.LengthFeasibleFlag = 0;
  }
  // 3. cost
  CostUpdate(Node_i);
}

void LddotUpperBound(const double &L_k, const double &Theta_k, const double &Thetadot_k, const double &Length_Low, const double &Length_Upp, double &Lddot_Low, double &Lddot_Upp)
{
  // This function is used to calculate the upper bound of the Lddot
  // The observation is that the maximum acceleration can only be achieved at the minimum height of the pendulum
  // Also this maximum acceleration is almost a constant in the upright direction, which is ~2.0g

  Lddot_Low = Clow * (L_k * Thetadot_k * Thetadot_k - g * cos(Theta_k));
  if(Length_Upp <=(Length_Low + L_unit))
  {
    Lddot_Upp = 0.0;
    return;
  }

  double Lddot_Upp_ref = Cupp * g * cos(Theta_k)/(Length_Low - Length_Upp) * (L_k - Length_Low) + Cupp * g * cos(Theta_k);
  Lddot_Upp = max(Lddot_Low, Lddot_Upp_ref);
  return;
}

std::set< std::tuple<int, int, int, int>> Backward_Propogation(const SystemIndex& Node_kp1_index, int & FeasiblePropFlag)
{
  // This function is used to conduct the backward propogation given a k+1 state
  // Here Length_Low/Upp is the actual pendulum length at the current setting so does g

  /*

      L_k = L_(k+1) - Ldot_(k+1) * delta_t
      Ldot_k = Ldot_(k+1) - L_(k+1) * theta_(k+1)^2 * delta_t + u_k
      theta_k = theta_(k+1) - thetadot_(k+1) * delta_t
      thetadot_k = (1 + 2 * Ldot_(k+1) * delta_t/L_(k+1)) * thetadot_(k+1) - g/L_(k+1) * sin(theta_(k+1)) & delta_t

  */

  double L_kp1=           L_vector[Node_kp1_index.L_index];
  double Ldot_kp1 =       Ldot_vector[Node_kp1_index.Ldot_index];
  double Theta_kp1 =      Theta_vector[Node_kp1_index.Theta_index];
  double Thetadot_kp1 =   Thetadot_vector[Node_kp1_index.Thetadot_index];

  // Discrete backward integration

  std::set< std::tuple<int, int, int, int>> Backward_Set;
  double L_k, Ldot_k, Theta_k, Thetadot_k;

  L_k = L_kp1 - Ldot_kp1 * delta_t;

  if((L_k<Length_Low)||(L_kp1<Length_Low))
  {
    FeasiblePropFlag = 0;
    return Backward_Set;
  }
  if((L_k>Length_Upp)||(L_kp1>Length_Upp))
  {
    FeasiblePropFlag = 0;
    return Backward_Set;
  }
  Theta_k = Theta_kp1 - Thetadot_kp1 * delta_t;
  Thetadot_k = (1 + 2 * Ldot_kp1 * delta_t/L_kp1) * Thetadot_kp1 - g/L_kp1 * sin(Theta_kp1) * delta_t;

  double Lddot_Low, Lddot_Upp;
  LddotUpperBound(L_kp1, Theta_kp1, Thetadot_kp1, Length_Low, Length_Upp, Lddot_Low, Lddot_Upp);
  std::vector<double> u_k_vec = linspace(Lddot_Low, Lddot_Upp, F_Grids);

  Ldot_k = Ldot_kp1 - u_k_vec[0] * delta_t;           // The only variant value

  SystemState x_k(L_k, Ldot_k, Theta_k, Thetadot_k);
  SystemIndex x_k_index = SystemState2StateIndex(x_k);

  Backward_Set.insert(make_tuple(x_k_index.L_index, x_k_index.Ldot_index, x_k_index.Theta_index, x_k_index.Thetadot_index));

  for (int i = 0; i < F_Grids-1; i++)
  {
    Ldot_k = Ldot_kp1 - u_k_vec[i+1] * delta_t;
    int Ldot_index = LdotValue2Index(Ldot_k);
    Backward_Set.insert(make_tuple(x_k_index.L_index, Ldot_index, x_k_index.Theta_index, x_k_index.Thetadot_index));
  }

  FeasiblePropFlag = 1;
  return Backward_Set;
}
SystemIndex FrontierPop(std::set< std::tuple<int, int, int, int>> &Frontier_Set)
{
  // This function is used to pop up the first element from the Frontier
  std::set<std::tuple<int, int, int, int>>::iterator Frontier_Set_First = Frontier_Set.begin();
  SystemIndex State_index(std::get<0>(*Frontier_Set_First),std::get<1>(*Frontier_Set_First),std::get<2>(*Frontier_Set_First),std::get<3>(*Frontier_Set_First));
  Frontier_Set.erase (Frontier_Set_First);
  return State_index;
}
void FrontierUpdate(const SystemState& Node_kp1_index, std::set< std::tuple<int, int, int, int>>& Frontier_Set, const std::set< std::tuple<int, int, int, int>>& Visited_Set)
{
   // This function is used to update the set element inside the Frontier.
   SystemIndex x_index = SystemState2Index(Node_kp1_index);

   int AddorNotFlag = Visited_Set.count(SystemIndex2Tuple(x_index));
   switch (AddorNotFlag) {
     case 0:
     Frontier_Set.insert(make_tuple(x_index.L_index, x_index.Ldot_index, x_index.Theta_index, x_index.Thetadot_index));
   }
}
void FrontierUpdate(const SystemIndex & Node_kp1_index, std::set< std::tuple<int, int, int, int>>& Backward_Set, std::set< std::tuple<int, int, int, int>>& Frontier_Set, Eigen::Tensor<DBNode,4>& StateNodeMatrix, std::vector<int>& NodeParentIndVector, std::vector<float>& NodeCostVector)
{
  // This function has two main tasks:
  // 1. Update the Frontier_Set with the Backward_Set information with respect to the Visited_Set
  // 2. Update the DataBase node in the StateNodeMatrix if there is an improvement of the objective function and Update
  DBNode* Node_kp1_ptr, *Node_k_ptr;
  Node_kp1_ptr = &StateNodeMatrix(Node_kp1_index.L_index, Node_kp1_index.Ldot_index, Node_kp1_index.Theta_index, Node_kp1_index.Thetadot_index);

  while(Backward_Set.size()>0)
  {
    SystemIndex Node_k_index = FrontierPop(Backward_Set);
    Node_k_ptr = &StateNodeMatrix(Node_k_index.L_index, Node_k_index.Ldot_index, Node_k_index.Theta_index, Node_k_index.Thetadot_index);
    switch (Node_k_ptr->ParentIndex)
    {
      case -2:
            // This means that this DataBase Node has not been instantiated before.
            // For Node_k's Update: ParentIndex and cost update according to the information of Parent node
            Node_k_ptr->ParentIndex = Node_kp1_ptr->SelfIndex;
            Node_k_ptr->cost = Node_kp1_ptr->cost;

            NodeParentIndVector[Node_k_ptr->SelfIndex] = Node_kp1_ptr->ParentIndex;
            NodeCostVector[Node_k_ptr->SelfIndex] = Node_kp1_ptr->cost;

            // For Node_kp1, we need to push_back this Node_k_ptr into its Children_Nodes
            Node_kp1_ptr->Children_Nodes.push_back(Node_k_ptr);

            // Then we add this Node_k into our Frontier_set.
            Frontier_Set.insert(make_tuple(Node_k_ptr->NodeIndex.L_index, Node_k_ptr->NodeIndex.Ldot_index, Node_k_ptr->NodeIndex.Theta_index, Node_k_ptr->NodeIndex.Thetadot_index));
      break;
    }
  }
}
void HJBDataWriter(const std::vector<int>& NodeParentIndVector, const std::vector<float>& NodeCostVector, const int & Angle_i)
{
  // // *. Node Transition
  // FILE * StateTransFile = NULL;
  // string StateTransFileName = "PVKNextInd" + Angle + ".bin";
  // StateTransFile = fopen(StateTransFileName, "wb");
  // fwrite(&NodeParentIndVector[0], sizeof(int), NodeParentIndVector.size(), StateTransFile);
  // fclose(StateTransFile);

  // *. Cost
  FILE * StateCostFile = NULL;
  string StateCostFileName = "PVKFailureMetric" + std::to_string(Angle_i) + ".bin";
  const char * StateCostFileName_Char = StateCostFileName.c_str();
  StateCostFile = fopen(StateCostFileName_Char, "wb");
  fwrite(&NodeCostVector[0], sizeof(float), NodeCostVector.size(), StateCostFile);
  fclose(StateCostFile);
}

SystemIndex VectorIndexToSystemIndex(const int & x)
{
  // This function is used to generate the SystemIndex from the vector index
  int a = Ldot_Grids * Theta_Grids * Thetadot_Grids;
  int b = Theta_Grids * Thetadot_Grids;
  int c = Thetadot_Grids;

  int i = floor(x/a);
  int j = floor((x - a * i)/b);
  int k = floor((x -a * i - b * j)/c);
  int l = x - a * i - b * j - c * k;

  SystemIndex x_systemindex(i,j,k,l);
  return x_systemindex;
}

void StateNodeMatrixInit(Eigen::Tensor<DBNode,4> &StateNodeMatrix, std::vector<float> &NodeCostVector)
{
  // This function is used to initialize the State Node Matrix
  DBNode * Node_ptr;
  for (int i = 0; i < L_Grids; i++)
  {
    for (int j = 0; j < Ldot_Grids; j++)
    {
      for (int k = 0; k < Theta_Grids; k++)
      {
        for (int l = 0; l < Thetadot_Grids; l++)
        {
          Node_ptr = &StateNodeMatrix(i,j,k,l);
          Node_ptr->InitUpdate(i, j, k, l, L_Grids, Ldot_Grids, Theta_Grids, Thetadot_Grids, L_vector, Ldot_vector, Theta_vector, Thetadot_vector);
          DBNodeTLcUpdate(*Node_ptr);
          NodeCostVector[Node_ptr->SelfIndex] = Node_ptr->cost;
        }
      }
    }
  }
}

void RankingVectorWriter(const std::vector<int>& NodeRankingVector)
{
  // *. Node Transition
  FILE * RankingVectorFile = NULL;
  RankingVectorFile = fopen("PVKRankingVector.bin", "wb");
  fwrite(&NodeRankingVector[0], sizeof(int), NodeRankingVector.size(), RankingVectorFile);
  fclose(RankingVectorFile);
}

std::vector<int> RankingVectorReader()
{
  FILE * RankingVectorFile = NULL;
  RankingVectorFile = fopen("PVKRankingVector.bin", "rb");
  std::vector<int> RankingVector(L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
  std::fread(&RankingVector[0], sizeof(int), RankingVector.size(), RankingVectorFile);
  fclose(RankingVectorFile);
  return RankingVector;
}

int main()
{
  std::vector<double> StateVectorSpecs;

  StateVectorSpecs.push_back(LLow);
  StateVectorSpecs.push_back(LUpp);
  StateVectorSpecs.push_back(LdotLow);
  StateVectorSpecs.push_back(LdotUpp);
  StateVectorSpecs.push_back(ThetaLow);
  StateVectorSpecs.push_back(ThetaUpp);
  StateVectorSpecs.push_back(ThetadotLow);
  StateVectorSpecs.push_back(ThetadotUpp);

  StateVectorSpecs.push_back(L_Grids * 1.0);
  StateVectorSpecs.push_back(Ldot_Grids * 1.0);
  StateVectorSpecs.push_back(Theta_Grids * 1.0);
  StateVectorSpecs.push_back(Thetadot_Grids * 1.0);

  const int AngleLow = 0;
  const int AngleUpp = 80;
  const int AngleDiff = 10;

  // This part needs to be updated with new information.
  StateVectorSpecs.push_back(AngleLow * 1.0);
  StateVectorSpecs.push_back(AngleUpp * 1.0);
  StateVectorSpecs.push_back(AngleDiff * 1.0);

  // Save the PVkDataSpecs first into PVkDataSpecs.bin
  FILE * StateVectorSpecsFile = NULL;
  StateVectorSpecsFile = fopen("PVkDataSpecs.bin", "wb");
  fwrite(&StateVectorSpecs[0], sizeof(double), StateVectorSpecs.size(), StateVectorSpecsFile);
  fclose(StateVectorSpecsFile);

  /*
    Main computation objects initialization
  */
  StateVectorSubs(L_vector, Ldot_vector, Theta_vector, Thetadot_vector);

  // However, there are still three dimensions to be discretized for completeness.
  // Also the gravitational vector needs to be discretized
  std::clock_t start;     double duration;      start = std::clock();

  for (int i = 0; i < (AngleUpp - AngleLow)/AngleDiff + 1; i++)
  {
    int Angle_i = AngleLow + AngleDiff * i;

    g = 9.81 * cos(Angle_i*1.0/180 * 3.1415926535897);
    cout<<g<<endl;

    // Here we are going one step forward towards completeness.
    // At the final time the velocities need to be enumerated as well.

    std::set< std::tuple<int, int, int, int>> Frontier_Set, Visited_Set, Infeasi_Set, Backward_Set;
    Eigen::Tensor<DBNode,4> StateNodeMatrix(L_Grids, Ldot_Grids, Theta_Grids, Thetadot_Grids);

    std::vector<int> NodeParentIndVector(L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
    std::vector<float> NodeCostVector(L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
    StateNodeMatrixInit(StateNodeMatrix, NodeCostVector);

    std::vector<int> RankingVectorObj = sort_indexes(NodeCostVector);
    RankingVectorWriter(RankingVectorObj);
    // std::vector<int> RankingVectorObj = RankingVectorReader();
    int FeasiblePropFlag;
    double CoveragePercentage = 0.0;
    int UsedPoint = 0;
    const int TotalNumber = L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids;
    int VisitedNumber = 0;
    // while(CoveragePercentage<0.99999999)     // Here we will run this whole algorithm until the coverage of the data base is within an acceptance range.
    while(VisitedNumber<TotalNumber)     // Here we will run this whole algorithm until the coverage of the data base is within an acceptance range.
    {
      SystemIndex NodeSystemIndex = VectorIndexToSystemIndex(RankingVectorObj[UsedPoint]);
      DBNode *Node_kp1_ptr = &StateNodeMatrix(NodeSystemIndex.L_index,NodeSystemIndex.Ldot_index,NodeSystemIndex.Theta_index,NodeSystemIndex.Thetadot_index);
      switch (Node_kp1_ptr->ParentIndex)
      {
        case -2:
        // ParentIndex Update
        Node_kp1_ptr->ParentIndex = -1;         // This means that the current node is not linked to any other parent node.
        NodeParentIndVector[Node_kp1_ptr->SelfIndex] = Node_kp1_ptr->ParentIndex;
        NodeCostVector[Node_kp1_ptr->SelfIndex] = Node_kp1_ptr->cost;
        // Backward Propogation from the current k+1 state
        switch (Node_kp1_ptr->LengthFeasibleFlag)
        {
          // In the case where this given point is feasible but its theta could be non-viable.
          case 1:
          {
            FrontierUpdate(Node_kp1_ptr->NodeState, Frontier_Set, Visited_Set);
            while(Frontier_Set.size()>0)
            {
              SystemIndex Node_kp1_index = FrontierPop(Frontier_Set);
              Backward_Set = Backward_Propogation(Node_kp1_index, FeasiblePropFlag);
              switch (FeasiblePropFlag)
              {
                case 1:
                {
                  // This means that the x_k_vec solution is feasible to be stacked into the Visited_Set
                  Visited_Set.insert(SystemIndex2Tuple(Node_kp1_index));
                  FrontierUpdate(Node_kp1_index, Backward_Set, Frontier_Set, StateNodeMatrix, NodeParentIndVector, NodeCostVector);
                  break;
                }
                default:
                // The value for the infeasible in the state matrix is given to be -1;
                Visited_Set.insert(SystemIndex2Tuple(Node_kp1_index));
                Infeasi_Set.insert(SystemIndex2Tuple(Node_kp1_index));
                break;
              }
            }
            break;
          }
        }
      }
      VisitedNumber = Visited_Set.size() + Infeasi_Set.size();
      // CoveragePercentage = 1.0 * VisitedNumber/(1.0 * L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
      // std::printf("CoveragePercentage: %f\n", CoveragePercentage);
      UsedPoint = UsedPoint + 1;
    }
    HJBDataWriter(NodeParentIndVector, NodeCostVector, Angle_i);
  }
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  printf ("Total running time: %f s\n", duration);
  return 0;
}
