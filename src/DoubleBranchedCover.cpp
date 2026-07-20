// Double branched cover invariant for planar knotoids.
// C++ port of double_branched_cover/two_fold_py3.py. See DoubleBranchedCover.hh.

#include "DoubleBranchedCover.hh"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::vector;

// Following PolynomialInvariant, only a timeout is signalled with a thrown
// std::runtime_error (interpreted as TIMEOUT by main.cpp); every other error
// prints a message and exits.
[[noreturn]] static void dbc_error(const std::string &msg)
{
  std::cerr << "*********************************************************" << std::endl;
  std::cerr << "ERROR " << msg << std::endl;
  std::cerr << "*********************************************************" << std::endl;
  std::exit(1);
}

namespace {

typedef std::array<int, 4> Crossing4;            // one PD crossing: 4 arc labels
typedef vector<Crossing4> PD;                    // planar diagram
typedef pair<int, int> Corner;                   // (label, sign)
typedef vector<vector<Corner> > RegionList;      // regions, each a list of corners

const int NONE = 0;   // sentinel for a "None" Gauss label (real labels are nonzero)

//----------------------------------------------------------------------------
// Build PD from the (open) Gauss code.
PD build_pd(const vector<int> &gauss, const vector<int> &signs)
{
  PD pd;
  set<int> checked;
  map<int, int> pos;
  for (int i = 0; i < (int)gauss.size(); i++) pos[gauss[i]] = i;
  for (int x = 0; x < (int)gauss.size(); x++)
    {
      int g = gauss[x];
      int inv = pos[-g];
      int a = std::abs(g);
      if (checked.count(a)) continue;
      int s = signs[a - 1];
      if (g > 0 && s == 1)       pd.push_back({inv, x + 1, inv + 1, x});
      else if (g > 0 && s == -1) pd.push_back({inv, x, inv + 1, x + 1});
      else if (g < 0 && s == 1)  pd.push_back({x, inv + 1, x + 1, inv});
      else if (g < 0 && s == -1) pd.push_back({x, inv, x + 1, inv + 1});
      checked.insert(a);
    }
  return pd;
}

//----------------------------------------------------------------------------
// Build PD of the cover (arc labels taken mod m = 2*max_label, cyclic).
PD build_cover_pd(const vector<int> &gauss, const vector<int> &signs)
{
  PD pd;
  set<int> checked;
  int max_n = *std::max_element(gauss.begin(), gauss.end());
  int m = 2 * max_n;
  map<int, int> pos;
  for (int i = 0; i < (int)gauss.size(); i++) pos[gauss[i]] = i;
  for (int x = 0; x < (int)gauss.size(); x++)
    {
      int g = gauss[x];
      int inv = pos[-g];
      int a = std::abs(g);
      if (checked.count(a)) continue;
      int s = signs[a - 1];
      if (g > 0 && s == 1)
        pd.push_back({inv % m, (x + 1) % m, (inv + 1) % m, x % m});
      else if (g > 0 && s == -1)
        pd.push_back({inv % m, x % m, (inv + 1) % m, (x + 1) % m});
      else if (g < 0 && s == 1)
        pd.push_back({x % m, (inv + 1) % m, (x + 1) % m, inv % m});
      else if (g < 0 && s == -1)
        pd.push_back({x % m, inv % m, (x + 1) % m, (inv + 1) % m});
      checked.insert(a);
    }
  return pd;
}

//----------------------------------------------------------------------------
// Trace the regions of the diagram.
RegionList regions(const PD &pd)
{
  set<Corner> checked;
  RegionList regionlist;
  int nb_endpoints = 0;
  map<int, vector<Corner> > occ;
  for (int j = 0; j < (int)pd.size(); j++)
    for (int l = 0; l < 4; l++)
      occ[pd[j][l]].push_back({j, l});
  for (int i = 0; i < (int)pd.size(); i++)
    for (int k = 0; k < 4; k++)
      {
        if (checked.count({i, k})) continue;
        vector<Corner> region;
        int i1 = i, k1 = k;
        while (true)
          {
            checked.insert({i1, k1});
            int d = pd[i1][k1] - pd[i1][(k1 + 2) % 4];
            if (d == 1 || d < -1) region.push_back({pd[i1][k1], 1});
            else if (d == -1 || d > 1) region.push_back({pd[i1][k1], -1});

            bool found = false;
            for (const Corner &jl : occ[pd[i1][k1]])
              if (!(i1 == jl.first && k1 == jl.second))
                { i1 = jl.first; k1 = (jl.second + 1) % 4; found = true; break; }
            if (!found)
              {
                nb_endpoints++;
                if (pd[i1][k1] > pd[i1][(k1 + 2) % 4]) region.push_back({pd[i1][k1], -1});
                else region.push_back({pd[i1][k1], 1});
                if (nb_endpoints > 2)
                  dbc_error("double_branched_cover: too many endpoints");
                k1 = (k1 + 1) % 4;
              }
            if (i1 == i && k1 == k) break;
          }
        regionlist.push_back(region);
      }
  return regionlist;
}

// max region index r with label a in region r AND label b in region r
int max_shared_region(const map<int, set<int> > &reg_of, int a, int b)
{
  int best = -1;
  map<int, set<int> >::const_iterator ita = reg_of.find(a), itb = reg_of.find(b);
  if (ita == reg_of.end() || itb == reg_of.end()) return best;
  for (int r : ita->second) if (itb->second.count(r)) best = std::max(best, r);
  return best;
}

//----------------------------------------------------------------------------
// Which region borders each crossing corner (by arc label).
vector<Crossing4> crossing_regions(const PD &pd, const RegionList &regions_list)
{
  map<int, set<int> > reg_of;
  for (int ri = 0; ri < (int)regions_list.size(); ri++)
    for (const Corner &c : regions_list[ri]) reg_of[c.first].insert(ri);
  vector<Crossing4> out;
  for (const Crossing4 &x : pd)
    out.push_back({max_shared_region(reg_of, x[0], x[1]),
                   max_shared_region(reg_of, x[1], x[2]),
                   max_shared_region(reg_of, x[2], x[3]),
                   max_shared_region(reg_of, x[3], x[0])});
  return out;
}

// Which region borders each crossing corner (by (label,sign), for the cover).
vector<Crossing4> crossing_regions_2(const PD &pd, const RegionList &regions_list)
{
  map<Corner, set<int> > reg_of;
  for (int ri = 0; ri < (int)regions_list.size(); ri++)
    for (const Corner &c : regions_list[ri]) reg_of[c].insert(ri);
  auto mx = [&](Corner a, Corner b) {
    int best = -1;
    map<Corner, set<int> >::const_iterator ita = reg_of.find(a), itb = reg_of.find(b);
    if (ita == reg_of.end() || itb == reg_of.end()) return best;
    for (int r : ita->second) if (itb->second.count(r)) best = std::max(best, r);
    return best;
  };
  vector<Crossing4> out;
  for (const Crossing4 &x : pd)
    {
      int r0, r1, r2, r3;
      int diff = x[1] - x[3];
      if (diff == 1 || diff < -1)
        {
          r0 = mx({x[0], 1}, {x[1], 1});   r1 = mx({x[1], -1}, {x[2], 1});
          r2 = mx({x[2], -1}, {x[3], -1}); r3 = mx({x[3], 1}, {x[0], -1});
        }
      else
        {
          r0 = mx({x[0], 1}, {x[1], -1});  r1 = mx({x[1], 1}, {x[2], 1});
          r2 = mx({x[2], -1}, {x[3], 1});  r3 = mx({x[3], -1}, {x[0], -1});
        }
      out.push_back({r0, r1, r2, r3});
    }
  return out;
}

//----------------------------------------------------------------------------
// A small undirected graph with integer edge weights (parallel edges collapsed,
// matching the networkx MultiGraph access pattern used in the reference).
struct Graph
{
  map<pair<int, int>, int> weight;   // canonical key (min,max) -> weight
  map<pair<int, int>, int> arc;      // optional arc label per edge
  set<int> nodes;

  static pair<int, int> key(int u, int v) { return u < v ? pair<int,int>(u, v) : pair<int,int>(v, u); }
  void add_node(int u) { nodes.insert(u); }
  bool has_edge(int u, int v) const { return weight.count(key(u, v)) > 0; }
  int  get_edge(int u, int v) const { return weight.at(key(u, v)); }
  void set_edge(int u, int v, int w, int a = -1)
  {
    weight[key(u, v)] = w;
    if (a != -1) arc[key(u, v)] = a;
    nodes.insert(u); nodes.insert(v);
  }
};

// Dijkstra shortest path (by summed weight); returns the node path.
vector<int> dijkstra_path(const Graph &g, int src, int dst)
{
  if (src == dst) return {src};
  map<int, int> adjbuild;
  map<int, vector<int> > adj;
  for (const auto &kv : g.weight)
    { adj[kv.first.first].push_back(kv.first.second);
      adj[kv.first.second].push_back(kv.first.first); }
  const int INF = 1 << 29;
  map<int, int> dist, prev;
  for (int n : g.nodes) dist[n] = INF;
  dist[src] = 0;
  std::priority_queue<pair<int, int>, vector<pair<int, int> >, std::greater<pair<int, int> > > pq;
  pq.push({0, src});
  while (!pq.empty())
    {
      pair<int, int> top = pq.top(); pq.pop();
      int d = top.first, u = top.second;
      if (d > dist[u]) continue;
      for (int v : adj[u])
        {
          int w = g.get_edge(u, v);
          if (dist[u] + w < dist[v])
            { dist[v] = dist[u] + w; prev[v] = u; pq.push({dist[v], v}); }
        }
    }
  if (dist.count(dst) == 0 || dist[dst] >= INF)
    dbc_error("double_branched_cover: region graph disconnected");
  vector<int> path;
  for (int at = dst; ; at = prev[at]) { path.push_back(at); if (at == src) break; }
  std::reverse(path.begin(), path.end());
  return path;
}

// summed edge weight along a shortest path (the "nesting" number)
int dijkstra_distance(const Graph &g, int src, int dst)
{
  vector<int> path = dijkstra_path(g, src, dst);
  int s = 0;
  for (int k = 0; k + 1 < (int)path.size(); k++) s += g.get_edge(path[k], path[k + 1]);
  return s;
}

//----------------------------------------------------------------------------
// The region-adjacency graph labelled by arcs (for finding cut arcs).
Graph build_arc_graph(const vector<Crossing4> &crs, const PD &pd)
{
  Graph g;
  int maxnode = 0;
  for (const Crossing4 &c : crs) for (int j = 0; j < 4; j++) maxnode = std::max(maxnode, c[j]);
  for (int n = 0; n < maxnode; n++) g.add_node(n);
  for (int i = 0; i < (int)pd.size(); i++)
    {
      const int *r = crs[i].data();
      const int arcs[4] = {pd[i][1], pd[i][2], pd[i][3], pd[i][0]};   // arc for pair (0,1),(1,2),(2,3),(3,0)
      const int pa[4] = {0, 1, 2, 3}, pb[4] = {1, 2, 3, 0};
      for (int e = 0; e < 4; e++)
        {
          int u = r[pa[e]], v = r[pb[e]];
          if (g.has_edge(u, v)) { /* weight is 1; reference never re-adds here */ }
          else if (u != v) g.set_edge(u, v, 1, arcs[e]);
        }
    }
  return g;
}

int find_outside(const vector<int> &outer, const RegionList &regions_list)
{
  set<int> target(outer.begin(), outer.end());
  for (int i = 0; i < (int)regions_list.size(); i++)
    {
      set<int> labs;
      for (const Corner &c : regions_list[i]) labs.insert(c.first);
      if (labs == target) return i;
    }
  return -1;
}

pair<int, int> find_endpoint_regions(const vector<int> &gauss, const RegionList &regions_list)
{
  int max_cross = *std::max_element(gauss.begin(), gauss.end());
  int head = -1, tail = -1;
  for (int i = 0; i < (int)regions_list.size(); i++)
    {
      set<int> labs;
      for (const Corner &c : regions_list[i]) labs.insert(c.first);
      if (labs.count(0)) head = i;
      if (labs.count(2 * max_cross)) tail = i;
    }
  return {head, tail};
}

// arcs crossed an odd number of times on paths from the outer region to the
// two endpoint regions
vector<int> find_cut(const Graph &g, int outside, pair<int, int> end_regs)
{
  map<int, int> cnt;
  for (int target : {end_regs.first, end_regs.second})
    {
      vector<int> path = dijkstra_path(g, outside, target);
      for (int k = 0; k + 1 < (int)path.size(); k++)
        cnt[g.arc.at(Graph::key(path[k], path[k + 1]))]++;
    }
  vector<int> cuts;
  for (const auto &kv : cnt) if (kv.second % 2 == 1) cuts.push_back(kv.first);
  return cuts;
}

//----------------------------------------------------------------------------
// Build the double branched cover Gauss code. Returns (labels, flags, signs).
struct Cover { vector<int> labels; vector<int> flags; vector<int> signs; };

Cover build_cover(const vector<int> &gauss, const vector<int> &signs_in, const vector<int> &cuts)
{
  vector<int> tempgc = gauss;
  vector<int> cs = cuts;
  std::sort(cs.rbegin(), cs.rend());
  for (int cut : cs) tempgc.insert(tempgc.begin() + cut, NONE);

  vector<int> l0, l1;                       // l1: 0/1 flag, -1 = None
  bool flag = false;
  for (int t : tempgc)
    {
      if (t == NONE) { l0.push_back(NONE); l1.push_back(-1); flag = !flag; }
      else           { l0.push_back(t);    l1.push_back(flag ? 1 : 0); }
    }
  // templift = lift + reversed(lift); flags in second half are flipped
  vector<int> t0 = l0, t1 = l1;
  for (int i = (int)l0.size() - 1; i >= 0; i--) t0.push_back(l0[i]);
  for (int i = (int)l1.size() - 1; i >= 0; i--)
    { int x = l1[i]; t1.push_back(x == -1 ? -1 : (x ? 0 : 1)); }

  auto find_idx = [&](int lab, int fl) {
    for (int i = 0; i < (int)t0.size(); i++) if (t0[i] == lab && t1[i] == fl) return i;
    return -1;
  };

  Cover cover;
  set<pair<int, int> > checked;
  for (int idx = 0; idx < (int)t0.size(); idx++)
    {
      int lab = t0[idx], fl = t1[idx];
      if (lab == NONE) continue;
      pair<int, int> ck(std::abs(lab), fl);
      if (checked.count(ck)) continue;
      int start = find_idx(lab, fl);
      int end = find_idx(-lab, fl);
      checked.insert(ck);
      int none_count = 0;
      for (int k = start + 1; k < end; k++) if (t0[k] == NONE) none_count++;
      int s = signs_in[std::abs(lab) - 1];
      cover.signs.push_back(none_count % 2 == 0 ? s : -s);
    }
  for (int i = 0; i < (int)t0.size(); i++) if (t0[i] != NONE) cover.labels.push_back(t0[i]);
  for (int i = 0; i < (int)t1.size(); i++) if (t1[i] != -1)   cover.flags.push_back(t1[i]);
  return cover;
}

// Map each original arc label to its lifted positions.
map<int, vector<int> > new_out_labels(const Cover &cover, int orig_gauss_len)
{
  vector<int> tmplft = cover.labels;
  tmplft.insert(tmplft.begin() + (int)tmplft.size() / 2, NONE);
  int denom = (int)tmplft.size() - 1;
  int mid = denom / 2;
  map<int, vector<int> > in_out;
  for (int i = 0; i <= orig_gauss_len; i++) in_out[i];      // keys 0..len
  set<int> arc_check;
  for (int val : tmplft)
    {
      if (arc_check.count(val)) continue;
      arc_check.insert(val);
      vector<int> indices;
      for (int i = 0; i < (int)tmplft.size(); i++) if (tmplft[i] == val) indices.push_back(i % denom);
      in_out[indices[0]] = indices;
      if (indices[0] == mid) in_out[indices[0]].push_back(indices[0]);
    }
  return in_out;
}

// Renumber cover labels into a sequential Gauss code. Returns (gauss, signs).
pair<vector<int>, vector<int> > renumber(const Cover &cover)
{
  vector<pair<int, int> > checked;
  vector<int> labels = cover.labels;
  for (int i = 0; i < (int)labels.size(); i++)
    {
      pair<int, int> ck(std::abs(labels[i]), cover.flags[i]);
      int idx = -1;
      for (int j = 0; j < (int)checked.size(); j++) if (checked[j] == ck) { idx = j; break; }
      if (idx == -1) { checked.push_back(ck); idx = (int)checked.size() - 1; }
      labels[i] = labels[i] < 0 ? -(idx + 1) : (idx + 1);
    }
  return {labels, cover.signs};
}

// Find the two cover regions whose boundaries together give the lifted outer arcs.
pair<vector<int>, vector<int> > find_in_out(const RegionList &cover_regions,
                                            const map<int, vector<int> > &new_out,
                                            const vector<int> &orig_outer)
{
  vector<int> out_in;
  for (int old : orig_outer)
    { const vector<int> &v = new_out.at(old); out_in.insert(out_in.end(), v.begin(), v.end()); }
  set<int> out_in_set(out_in.begin(), out_in.end());

  RegionList all_regions;
  for (const vector<Corner> &region : cover_regions)
    {
      bool inside = true;
      for (const Corner &c : region) if (!out_in_set.count(c.first)) { inside = false; break; }
      if (inside) all_regions.push_back(region);
    }
  RegionList both;
  int flag = 0;
  for (int k = 0; k < (int)all_regions.size(); k++)
    for (int l = k + 1; l < (int)all_regions.size(); l++)
      {
        set<int> u;
        for (const Corner &c : all_regions[k]) u.insert(c.first);
        for (const Corner &c : all_regions[l]) u.insert(c.first);
        if (u == out_in_set) { both.push_back(all_regions[k]); both.push_back(all_regions[l]); flag++; }
      }
  if (flag > 1) dbc_error("double_branched_cover: more than two boundary regions");
  if (both.size() < 2) dbc_error("double_branched_cover: boundary regions not found");
  vector<int> outside, inside_;
  for (const Corner &c : both[0]) outside.push_back(c.first);
  for (const Corner &c : both[1]) inside_.push_back(c.first);
  return {outside, inside_};
}

//----------------------------------------------------------------------------
// State connectivity graph for one Kauffman state (weights 0/1/2).
Graph build_graph_2(const vector<int> &state_bits, const vector<Crossing4> &crs)
{
  Graph g;
  int node_count = 0;
  for (const Crossing4 &c : crs) for (int j = 0; j < 4; j++) node_count = std::max(node_count, c[j]);
  for (int n = 0; n < node_count; n++) g.add_node(n);

  auto corner = [&](int u, int v) {
    if (g.has_edge(u, v)) { if (g.get_edge(u, v) > 1) g.set_edge(u, v, 1); }
    else if (u != v) g.set_edge(u, v, 1);
  };
  for (int i = 0; i < (int)state_bits.size(); i++)
    {
      int c0 = crs[i][0], c1 = crs[i][1], c2 = crs[i][2], c3 = crs[i][3];
      corner(c0, c1); corner(c1, c2); corner(c2, c3); corner(c3, c0);
      if (state_bits[i] == 0)
        {
          if (!g.has_edge(c0, c2) && c0 != c2) g.set_edge(c0, c2, 2);
          if (g.has_edge(c1, c3)) { if (g.get_edge(c1, c3) > 0) g.set_edge(c1, c3, 0); }
          else if (c1 != c3) g.set_edge(c1, c3, 0);
        }
      else
        {
          if (g.has_edge(c0, c2)) { if (g.get_edge(c0, c2) > 0) g.set_edge(c0, c2, 0); }
          else if (c0 != c2) g.set_edge(c0, c2, 0);
          if (!g.has_edge(c1, c3) && c1 != c3) g.set_edge(c1, c3, 2);
        }
    }
  return g;
}

// number of connected components of the weight-0 subgraph, minus one
int component_counter(const Graph &g)
{
  map<int, int> parent;
  for (int n : g.nodes) parent[n] = n;
  std::function<int(int)> findp = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
  for (const auto &kv : g.weight)
    if (kv.second == 0)
      { int a = findp(kv.first.first), b = findp(kv.first.second); if (a != b) parent[a] = b; }
  set<int> roots;
  for (int n : g.nodes) roots.insert(findp(n));
  return (int)roots.size() - 1;
}

int writhe_2(const PD &pd)
{
  int writhe = 0;
  for (const Crossing4 &x : pd)
    {
      int d = x[1] - x[3];
      if (d == 1 || d < -1) writhe += 1;
      else if (d == -1 || d > 1) writhe -= 1;
    }
  return writhe;
}

int find_region_by_labels(const vector<int> &labels, const RegionList &regions_list)
{
  set<int> target(labels.begin(), labels.end());
  for (int i = 0; i < (int)regions_list.size(); i++)
    {
      set<int> labs;
      for (const Corner &c : regions_list[i]) labs.insert(c.first);
      if (labs == target) return i;
    }
  return -1;
}

} // anonymous namespace

//----------------------------------------------------------------------------
// Public: build the double branched cover (gc_in_st).
DBCLift double_branched_cover_lift(const vector<int> &gauss,
                                   const vector<int> &signs,
                                   const vector<int> &outer)
{
  PD pd = build_pd(gauss, signs);
  RegionList regions_list = regions(pd);
  int outside_region = find_outside(outer, regions_list);
  if (outside_region < 0)
    dbc_error("double_branched_cover: outer arcs do not match any region");
  pair<int, int> end_regs = find_endpoint_regions(gauss, regions_list);
  vector<Crossing4> crs_reg = crossing_regions(pd, regions_list);
  Graph arc_graph = build_arc_graph(crs_reg, pd);
  vector<int> cuts = find_cut(arc_graph, outside_region, end_regs);
  Cover cover = build_cover(gauss, signs, cuts);
  map<int, vector<int> > new_out = new_out_labels(cover, (int)gauss.size());
  pair<vector<int>, vector<int> > covergc = renumber(cover);
  PD coverpd = build_cover_pd(covergc.first, covergc.second);
  RegionList cover_regions = regions(coverpd);
  pair<vector<int>, vector<int> > in_out = find_in_out(cover_regions, new_out, outer);

  DBCLift lift;
  lift.gauss = covergc.first;
  lift.signs = covergc.second;
  lift.outside = in_out.first;
  lift.inside = in_out.second;
  return lift;
}

//----------------------------------------------------------------------------
// Public: Jones polynomial of a knot in a solid torus (jones_poly_st).
Polynomial solid_torus_jones_polynomial(const DBCLift &lift, time_t timeout)
{
  const vector<std::string> vars = {"A", "v"};
  PD pd = build_cover_pd(lift.gauss, lift.signs);
  RegionList regions_list = regions(pd);
  int writhe = writhe_2(pd);
  int outside_region = find_region_by_labels(lift.outside, regions_list);
  int inside_region = find_region_by_labels(lift.inside, regions_list);
  if (outside_region < 0 || inside_region < 0)
    dbc_error("double_branched_cover: lift boundary region not found");
  vector<Crossing4> crs_reg = crossing_regions_2(pd, regions_list);

  int n = (int)pd.size();

  // loop factor (-A^2 - A^-2)
  Polynomial loop(vars);
  loop.add(-1.0, "A", 2);
  loop.add(-1.0, "A", -2);

  Polynomial total(vars);
  time_t t0 = time(NULL);
  for (long mask = 0; mask < (1L << n); mask++)
    {
      if (timeout > 0 && time(NULL) - t0 >= timeout)
        throw std::runtime_error("DoubleBranchedCover: Timeout.");
      vector<int> bits(n);
      int ones = 0;
      for (int i = 0; i < n; i++) { bits[i] = (mask >> i) & 1; ones += bits[i]; }
      int a_exp = (n - ones) - ones;   // A^(#0 - #1)

      Graph g = build_graph_2(bits, crs_reg);
      int comp = component_counter(g);
      int nest = (outside_region == inside_region) ? 0
                 : dijkstra_distance(g, outside_region, inside_region);
      int delta = (nest != 0) ? (comp - nest) : (comp - 1);
      if (delta < 0)
        dbc_error("double_branched_cover: negative loop exponent");

      // term = A^a_exp * loop^delta * v^nest
      Polynomial term(vars);
      term.add(1.0, "A", a_exp);
      for (int d = 0; d < delta; d++) term.multiply(loop);
      if (nest != 0) { Polynomial vpow(vars); vpow.add(1.0, "v", nest); term.multiply(vpow); }
      total.add(term);
    }

  // multiply by (-A^3)^(-writhe) = (-1)^writhe * A^(-3*writhe)
  Polynomial norm(vars);
  norm.add(writhe % 2 == 0 ? 1.0 : -1.0, "A", -3 * writhe);
  total.multiply(norm);
  return total;
}

//----------------------------------------------------------------------------
Polynomial double_branched_cover_polynomial(const vector<int> &gauss,
                                            const vector<int> &signs,
                                            const vector<int> &outer,
                                            time_t timeout)
{
  return solid_torus_jones_polynomial(double_branched_cover_lift(gauss, signs, outer), timeout);
}

//----------------------------------------------------------------------------
// Public: full invariant from a KnotoID planar PlanarDiagram. Extracts the
// extended oriented Gauss code (crossings + signs + exterior arcs) produced by
// PlanarDiagram::save_to_file_extended_gauss_code and feeds it to the raw
// Gauss-code implementation. Arcs are numbered in traversal order by that
// writer, matching the position-based arc labels this invariant expects.
Polynomial double_branched_cover_polynomial(PlanarDiagram &diagram, time_t timeout)
{
  std::ostringstream oss;
  diagram.save_to_file_extended_gauss_code(oss, true);   // one line

  vector<int> gauss, signs, outer;
  bool signs_seen = false;
  std::istringstream iss(oss.str());
  std::string tok;
  while (iss >> tok)
    {
      bool has_digit = tok.find_first_of("0123456789") != std::string::npos;
      if (!has_digit)                       // token of only +/- => the sign string
        {
          for (char c : tok) signs.push_back(c == '+' ? 1 : -1);
          signs_seen = true;
        }
      else if (!signs_seen)
        gauss.push_back(std::stoi(tok));    // signed crossing sequence
      else
        outer.push_back(std::stoi(tok));    // exterior arc labels
    }

  if (gauss.empty())
    dbc_error("double_branched_cover: empty or non-knotoid diagram");
  if (signs.empty())
    dbc_error("double_branched_cover: could not read crossing signs");
  return double_branched_cover_polynomial(gauss, signs, outer, timeout);
}
