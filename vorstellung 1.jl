using LinearAlgebra
using Pkg
using Plots

import Plots: plot, plot!, scatter, scatter!
Plots.gr()  # Initialisiert das GR-Backend für die Plot-Bibliothek

struct Point
    x::Float64
    y::Float64
end
# Abstrakter Typ für ein Face (z. B. Dreieck), von dem weitere Strukturen erben können
abstract type Face end

# Kantenstruktur in der DCEL-Datenstruktur
mutable struct Edge
    origin::Point
    twin::Union{Edge, Nothing}
    next::Union{Edge, Nothing}
    prev::Union{Edge, Nothing}
    face::Union{Face, Nothing}
end

# Repräsentiert ein Dreieck, das aus einer seiner Kanten besteht
mutable struct Dreieck <: Face
    edge::Union{Edge, Nothing}
end

# Struktur zur Speicherung der Delaunay-Triangulation
mutable struct Delaunay
    triangles::Set{Dreieck}
    points::Set{Point}
    bounding_triangle::Union{Dreieck, Nothing}
    Delaunay() = new(Set{Dreieck}(), Set{Point}(), nothing)
end

# Repräsentiert das Ergebnis eines Voronoi-Diagramms
# Kanten des Voronoi-Diagramms
 # Zuordnung von Punkten zu ihren Regionen
struct VoronoiDiagram
    edges::Vector{Tuple{Point, Point}}
    regions::Dict{Point, Vector{Point}}
end

 # Setzt zyklische Nachbarschaftsbeziehungen der Kanten im Dreieck
 # Erstellt ein Dreieck mit einer beliebigen seiner Kanten
 # Verknüpft jede Kante mit dem zugehörigen Dreieck (Face-Zuordnung)
function create_triangle(a, b, c, e1, e2, e3)
    e1.next = e2; e1.prev = e3
    e2.next = e3; e2.prev = e1
    e3.next = e1; e3.prev = e2
    tri = Dreieck(e1)
    e1.face = tri; e2.face = tri; e3.face = tri
    return tri
end

#Diese Funktion prüft, ob der gegenüberliegende Punkt d außerhalb oder auf dem Umkreis des Dreiecks (a,b,c) liegt.
# Punkte des aktuellen Dreiecks
# Wenn es keine gegenüberliegende Kante gibt (Randkante),
# gilt die Kante als gültig (kein Flip notwendig)
function check_umkreis(e::Edge)::Bool
    a = e.origin
    b = e.next.origin
    c = e.next.next.origin

    if e.twin === nothing || e.twin.face === nothing
        return true
    end

    d = e.twin.next.next.origin

    M = [a.x a.y a.x^2 + a.y^2 1;
         b.x b.y b.x^2 + b.y^2 1;
         c.x c.y c.x^2 + c.y^2 1;
         d.x d.y d.x^2 + d.y^2 1]

    return det(M) <= 0 #nicht in Umkreis and auf Umkereis
end


#Diese Funktion „flip！“ eine Kante in der Delaunay-Triangulation, falls sie nicht Delaunay-konform ist.
# Prüfen, ob die Kante oder ihre Zwillinge existieren
function flip!(e::Edge, D::Delaunay)
    # Basic check: ensure the edge and its twin exist and have faces
    # If any are missing, a flip can't happen, so we return early.
    if e.twin === nothing || e.face === nothing || e.twin.face === nothing
        return D, e # Return the Delaunay structure and the original edge
    end

    #  Identify the four points of the quadrilateral (a-c-b-d)
    # The edge 'e' goes from 'a' to 'c' within tri1.
    # The edge 'e.twin' goes from 'c' to 'a' within tri2.
    a = e.origin
    c = e.next.origin # Point 'c' is the target of edge 'e'
    b = e.next.next.origin # Point 'b' is the remaining point in tri1 (after a, c)
                           # e.prev.origin (origin of edge b->a)
    d = e.twin.next.next.origin # Point 'd' is the remaining point in tri2 (after c, a)
                                # e.twin.prev.origin (origin of edge d->c)

    # Get the two old triangles involved in the flip
    tri1_old = e.face      # Triangle (a, c, b)
    tri2_old = e.twin.face # Triangle (c, a, d)

    # Identify the four EXISTING edges that form the outer boundary of the quadrilateral (a-b, b-c, c-d, d-a)
    # These edges will be reused. We just need to update their internal pointers and faces.

    # Edges from tri1_old:
    edge_ac_old = e              # The edge being flipped (a -> c)
    edge_cb_old = e.next         # The edge (c -> b)
    edge_ba_old = e.prev         # The edge (b -> a)

    # Edges from tri2_old:
    edge_ca_old = e.twin         # The twin of the edge being flipped (c -> a)
    edge_ad_old = e.twin.next    # The edge (a -> d)
    edge_dc_old = e.twin.prev    # The edge (d -> c)


    #  Remove the old triangles from the Delaunay structure's set of triangles
    delete!(D.triangles, tri1_old)
    delete!(D.triangles, tri2_old)

    # Create the NEW shared edge (b, d) and its twin (d, b)
    # This is the ONLY place where new Edge objects are created in flip!
    new_e_bd = Edge(b, nothing, nothing, nothing, nothing) # Edge from B to D
    new_e_db = Edge(d, nothing, nothing, nothing, nothing) # Edge from D to B (twin of new_e_bd)

    new_e_bd.twin = new_e_db
    new_e_db.twin = new_e_bd


    tri1_new = create_triangle(a, d, b, edge_ad_old, new_e_db, edge_ba_old)

    # New Tri 2: (d, c, b)
    # Edges: (d->c) (existing: edge_dc_old), (c->b) (existing: edge_cb_old), (b->d) (new: new_e_bd)
    tri2_new = create_triangle(d, c, b, edge_dc_old, edge_cb_old, new_e_bd)

    #  Update the `face` pointers for the new central edges
    new_e_bd.face = tri2_new # Edge b->d belongs to tri2_new
    new_e_db.face = tri1_new # Edge d->b belongs to tri1_new



    if edge_ba_old.twin !== nothing
        edge_ba_old.twin.twin = edge_ba_old
    end
    if edge_cb_old.twin !== nothing
        edge_cb_old.twin.twin = edge_cb_old
    end
    if edge_ad_old.twin !== nothing
        edge_ad_old.twin.twin = edge_ad_old
    end
    if edge_dc_old.twin !== nothing
        edge_dc_old.twin.twin = edge_dc_old
    end

    # Add the new triangles to the Delaunay structure
    push!(D.triangles, tri1_new)
    push!(D.triangles, tri2_new)

    # Return the Delaunay structure and the new central edge (new_e_bd)
    # This new edge (b->d) is the one that might need further recursive flips.
    return D, new_e_bd
end

#Dieser Algorithmus prüft, ob eine Kante die Delaunay-Bedingung verletzt (also ob ein gegenüberliegender Punkt im Umkreis eines Dreiecks liegt). Falls ja, wird die Kante geflippt, und der Algorithmus ruft sich rekursiv auf den neuen angrenzenden Kanten auf. So wird Schritt für Schritt das gesamte Netz Delaunay-konform gemacht.
function recursive_flip!(e::Edge, D::Delaunay)
    if !check_umkreis(e)
        D, new_e = flip!(e, D)

        if new_e.next !== nothing && new_e.next.twin !== nothing
            recursive_flip!(new_e.next.twin, D)
        end
        if new_e.prev !== nothing && new_e.prev.twin !== nothing
            recursive_flip!(new_e.prev.twin, D)
        end
    end #prüfen jede Edge von e.face
end

#Hilfsfunktion, Prüft, ob ein Punkt p innerhalb eines gegebenen Dreiecks tri liegt
function point_in_triangle(p::Point, tri::Dreieck)
    a = tri.edge.origin
    b = tri.edge.next.origin
    c = tri.edge.next.next.origin
    # Berechne Richtungsvektoren des Dreiecks
    v0 = [c.x - a.x, c.y - a.y]
    v1 = [b.x - a.x, b.y - a.y]
    v2 = [p.x - a.x, p.y - a.y]

 # Berechne den Nenner der baryzentrischen Koordinatenformel
 # Falls der Nenner fast null ist, ist das Dreieck entartet
    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1)
    dot02 = dot(v0, v2)
    dot11 = dot(v1, v1)
    dot12 = dot(v1, v2)

    denom = dot00 * dot11 - dot01^2
    if abs(denom) < 1e-10
        return false
    end

     # Berechne die baryzentrischen Koordinaten (u, v)
    inv_denom = 1 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom

    return u >= 0 && v >= 0 && u + v <= 1
end

#Hilfsfunktion, # Findet das Dreieck in der Delaunay-Triangulation D, das den Punkt p enthält
function find_containing_triangle(p::Point, D::Delaunay)
    for tri in D.triangles
        point_in_triangle(p, tri) && return tri
    end
    return nothing
end

# Füge den neuen Punkt zur Punktmenge hinzu
function insert_point!(p::Point, D::Delaunay)
    push!(D.points, p)
    # Füge den neuen Punkt zur Punktmenge hinzu
    #bilden eine bounding_triangle,wie oben
    if isempty(D.triangles)
        n = 1000.0
        a = Point(-n, -n)
        b = Point(3n, -n)
        c = Point(0, 3n)

        e1 = Edge(a, nothing, nothing, nothing, nothing)
        e2 = Edge(b, nothing, nothing, nothing, nothing)
        e3 = Edge(c, nothing, nothing, nothing, nothing)
        
        tri = create_triangle(a, b, c, e1, e2, e3)
        D.bounding_triangle = tri
        push!(D.triangles, tri)
    end

    # Finde das Dreieck, das den Punkt p enthält
    T = find_containing_triangle(p, D)
    if T === nothing
        println("Der Punkt ($(p.x), $(p.y)) liegt in keinem Dreieck.")
        return D
    end
    
     # Kanten des Dreiecks T
    e1 = T.edge
    e2 = e1.next
    e3 = e2.next
    a, b, c = e1.origin, e2.origin, e3.origin

    e_ap = Edge(a, nothing, nothing, nothing, nothing)
    e_pa = Edge(p, nothing, nothing, nothing, nothing)
    e_bp = Edge(b, nothing, nothing, nothing, nothing)
    e_pb = Edge(p, nothing, nothing, nothing, nothing)
    e_cp = Edge(c, nothing, nothing, nothing, nothing)
    e_pc = Edge(p, nothing, nothing, nothing, nothing)

    e_ap.twin = e_pa; e_pa.twin = e_ap
    e_bp.twin = e_pb; e_pb.twin = e_bp
    e_cp.twin = e_pc; e_pc.twin = e_cp

    tri1 = create_triangle(a, b, p, e1, e_bp, e_pa)
    tri2 = create_triangle(b, c, p, e2, e_cp, e_pb)
    tri3 = create_triangle(c, a, p, e3, e_ap, e_pc)
    
    if e1.twin !== nothing; e1.twin.twin = e1; end
    if e2.twin !== nothing; e2.twin.twin = e2; end
    if e3.twin !== nothing; e3.twin.twin = e3; end

    delete!(D.triangles, T)
    push!(D.triangles, tri1, tri2, tri3)

    if e1.twin !== nothing; recursive_flip!(e1.twin, D); end
    if e2.twin !== nothing; recursive_flip!(e2.twin, D); end
    if e3.twin !== nothing; recursive_flip!(e3.twin, D); end

    return D
end

# Überprüft, ob zwei Punkte p1 und p2 gleich sind (innerhalb einer Toleranzgrenze).
# Da bei Gleitkomma-Berechnungen kleine Rundungsfehler auftreten können,
# wird statt eines exakten Vergleichs eine Toleranz verwendet.
function point_equals(p1::Point, p2::Point, tolerance=1e-6)
    return abs(p1.x - p2.x) < tolerance && abs(p1.y - p2.y) < tolerance
end

#Voronoi diagram Brechnen
#Voronoi diagram Brechnen
function compute_voronoi(D::Delaunay)
    voronoi_edges = Vector{Tuple{Point, Point}}()
    circumcenters = Dict{Dreieck, Point}()
    regions = Dict{Point, Vector{Point}}()

    #Umkreismittelpunkte berechnen und Regionen aufbauen
    for tri in D.triangles
        # *** FIRST CRITICAL ADDITION: Skip the bounding triangle for circumcenter calculation ***
        # Its circumcenter is usually very far away and not part of the actual Voronoi diagram of user points.
        if tri === D.bounding_triangle
            continue
        end

        edge = tri.edge
        a = edge.origin
        b = edge.next.origin
        c = edge.next.next.origin

        d = 2 * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y))
        if abs(d) < 1e-9
            continue # Überspringe entartete Dreiecke (collinear points)
        end

        # Berechnung der Koordinaten des Umkreismittelpunktes
        ux = (a.x^2 + a.y^2) * (b.y - c.y) +
             (b.x^2 + b.y^2) * (c.y - a.y) +
             (c.x^2 + c.y^2) * (a.y - b.y)
        uy = (a.x^2 + a.y^2) * (c.x - b.x) +
             (b.x^2 + b.y^2) * (a.x - c.x) +
             (c.x^2 + c.y^2) * (b.x - a.x)

        cc = Point(ux / d, uy / d) # Umkreismittelpunkt
        circumcenters[tri] = cc

        # Füge diesen Mittelpunkt zu jeder zugehörigen Punktregion hinzu
        # (This part primarily affects region drawing, less so line messiness)
        for p in (a, b, c)
            # You might also want to skip bounding triangle points here for region building.
            # However, for messy lines, the main fix is in the circumcenter and edge loop below.
            if !haskey(regions, p)
                regions[p] = Point[]
            end
            push!(regions[p], cc)
        end
    end

    # Voronoi Kanten basierend auf benachbarten Umkreismittelpunkten erzeugen
    for tri in D.triangles
        # *** SECOND CRITICAL ADDITION: Skip the bounding triangle itself from processing for edges ***
        # We don't want to draw edges *from* its circumcenter.
        if tri === D.bounding_triangle
            continue
        end

        if !haskey(circumcenters, tri)
            continue # Should not happen if bounding triangle skipped above
        end
        cc1 = circumcenters[tri]

        edge = tri.edge
        for _ in 1:3
            t = edge.twin
            # Ensure twin edge exists, its face exists, and its face's circumcenter is calculated.
            if t !== nothing && t.face !== nothing && haskey(circumcenters, t.face)
                # *** THIRD CRITICAL ADDITION: If the twin triangle is the bounding triangle,
                # this Voronoi edge extends to infinity. We do NOT connect to its
                # (potentially distant) circumcenter here.
                if t.face === D.bounding_triangle
                    edge = edge.next
                    continue # Skip this edge, its Voronoi dual is infinite
                end

                cc2 = circumcenters[t.face]

                # um doppelte Kanten zu vermeiden (lexicographical order)
                if cc1.x < cc2.x || (abs(cc1.x - cc2.x) < 1e-9 && cc1.y < cc2.y)
                    push!(voronoi_edges, (cc1, cc2))
                end
            end
            edge = edge.next
        end
    end

    return VoronoiDiagram(voronoi_edges, regions)
end


# ==================GAME========================

struct GamePoint
    point::Point
    owner::Int
    radius::Float64 #für Zeichnen
end

# Konstruktor für GamePoint: Erzeugt einen Punkt mit gegebenen Koordinaten und Spieler
GamePoint(x::Real, y::Real, owner::Int) = GamePoint(Point(Float64(x), Float64(y)), owner, 6.0)

# Initialisierung der Spielzustände
game_points = GamePoint[]  # Liste aller gesetzten Spielpunkte
current_player = 1  # Aktuell aktiver Spieler (Spieler 1 beginnt)
turn_count = 0   # Anzahl der bisherigen Züge
max_turns = 20    # Maximale Anzahl der Spielzüge
game_over = false  # Spielstatus: false = läuft noch
winner = nothing   # Gewinner (noch nicht bestimmt)
canvas_size = 100.0  # Größe der Zeichenfläche

# Initialisierung der Delaunay-Triangulation (für die Punktverteilung)
delaunay_triangulation = Delaunay()
current_voronoi = nothing # Voronoi-Diagramm

# ===============Funktion von game ===========================

#Funktion zum Zurücksetzen des Spiels
function reset_game()
    global game_points, current_player, turn_count, game_over, winner
    global delaunay_triangulation, current_voronoi
    
    empty!(game_points)  # Alle Spielpunkte löschen
    current_player = 1  # Setze aktuellen Spieler auf Spieler 1 zurück
    turn_count = 0    # Setze Zugzähler zurück
    game_over = false  # Spielstatus auf „nicht beendet“ setzen
    winner = nothing  # Gewinner zurücksetzen...wie oben
    delaunay_triangulation = Delaunay()
    current_voronoi = nothing
    
    println("Das Spiel wurde zurückgesetzt.")
end


function add_point(x::Real, y::Real)
    global game_points, current_player, turn_count, game_over, winner
    global delaunay_triangulation, current_voronoi
    
    # Eingabekoordinaten in Float64 umwandeln
    x_f, y_f = Float64(x), Float64(y)
    
    # Prüfen, ob der Punkt innerhalb der Spielfeldgrenzen liegt
    if x_f < 0 || x_f > canvas_size || y_f < 0 || y_f > canvas_size
        println("Der angeklickte Ort liegt außerhalb der Grenzen.")
        return false
    end
    
     # Versuche, einen nahegelegenen Punkt zu entfernen (z.B. zum Überschreiben)
    if remove_nearby_point(x_f, y_f)
        rebuild_delaunay_and_voronoi()  # Aktualisiere Delaunay-Triangulation und Voronoi-Diagramm
        if game_over
            winner = determine_winner() # Bestimme den Gewinner, falls das Spiel beendet ist
        end
        show_game() # Zeige den aktuellen Spielstand
        return true
    end
    
    if game_over
        println("Das Spiel ist beendet.Bitte führen Sie restart() aus, um es neu zu beginnen.")
        return false
    end
    
    new_game_point = GamePoint(x_f, y_f, current_player)
    push!(game_points, new_game_point)
    
    insert_point!(new_game_point.point, delaunay_triangulation)
    current_voronoi = compute_voronoi(delaunay_triangulation)
    
    println("Player $current_player in ($x, $y) platziert einen Punkt")
    
    turn_count += 1
    
    if turn_count >= max_turns
        game_over = true
        winner = determine_winner() #unter
        if winner === nothing
            println("Spiel beendet - Unentschieden!")
        else
            println("Spiel beendet - Player $winner siegt.")
        end
    else
        # Spieler wechseln
        current_player = (current_player == 1) ? 2 : 1
        println("jetzt $current_player")
    end
    
    show_game()
    return true
end

#Diese Funktion remove_nearby_point überprüft, ob sich in der Nähe der angegebenen Koordinaten (x, y) ein Spielpunkt (GamePoint) befindet. Falls ja, wird dieser Punkt gelöscht.
function remove_nearby_point(x::Float64, y::Float64)
    global game_points
    click_threshold_sq = 0.0

    # Durchlaufe alle Spielpunkte
    for (i, gp) in enumerate(game_points)
        # Schwellenwert als Quadrat des Radius des Punkts
        click_threshold_sq = gp.radius^2 #unter
        if (x - gp.point.x)^2 + (y - gp.point.y)^2 < click_threshold_sq
            deleted_gp = popat!(game_points, i) #delete point
            println("Punkt von Spieler  $(deleted_gp.owner) wurde gelöscht.")
            return true
        end
    end
    
    return false
end

function rebuild_delaunay_and_voronoi()
    global delaunay_triangulation, current_voronoi, game_points
    
    delaunay_triangulation = Delaunay()
     # Füge alle Punkte aus game_points in die Delaunay-Triangulation ein
    for gp in game_points
        insert_point!(gp.point, delaunay_triangulation)
    end

    # Wenn Punkte vorhanden sind, berechne das Voronoi-Diagramm
    if !isempty(game_points)
        current_voronoi = compute_voronoi(delaunay_triangulation)
    else
        # Falls keine Punkte vorhanden sind, setze current_voronoi auf nichts
        current_voronoi = nothing
    end
end

function find_closest_game_point(test_point::Point) #mit remove_nearby_point normalerweise, verminden
    global game_points
    
    if isempty(game_points)
        return nothing
    end
    
    min_dist_sq = Inf
    closest_gp = nothing
    
    for gp in game_points
        dist_sq = (test_point.x - gp.point.x)^2 + (test_point.y - gp.point.y)^2
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq
            closest_gp = gp
        end
    end
    
    return closest_gp
end

function calculate_areas_voronoi()
    step = 1.0 # Rastergröße für die Flächenabschätzung
    area1 = 0.0
    area2 = 0.0

     # Durchlaufe alle Punkte im Raster auf der Zeichenfläche
    for x in 0:step:canvas_size
        for y in 0:step:canvas_size
            test_point = Point(x, y)
            closest_game_point = find_closest_game_point(test_point)
            
            if closest_game_point !== nothing
                # Fläche je nach Eigentümer zuweisen
                if closest_game_point.owner == 1
                    area1 += step * step
                else
                    area2 += step * step
                end
            end
        end
    end
    
    return area1, area2
end

function determine_winner()
    area1, area2 = calculate_areas_voronoi()
    
    if abs(area1 - area2) < 1.0
        return nothing
    elseif area1 > area2
        return 1
    else
        return 2
    end
end

# ====================Zeichnen======================

function show_game()
    global game_points, current_voronoi, delaunay_triangulation
    
    p = Plots.plot(
        xlim=(0, canvas_size),
        ylim=(0, canvas_size),
        aspect_ratio=:equal,# gleiche Skalierung auf beiden Achsen
        title=get_game_title(),
        size=(700, 700),
        legend=:topright,
        grid=false,
        ticks=false,
        border=:none  # Kein Rand um den Plot
    )
    
    # Falls Spielpunkte vorhanden sind, zeichne die Voronoi-Regionen (pixelbasiert)
    if !isempty(game_points)
        draw_voronoi_regions_pixel_based!(p)
    end
    

    if current_voronoi !== nothing
        draw_voronoi_boundaries!(p, current_voronoi)
    end
    
    draw_player_points!(p)
    display(p)
    return p
end

function draw_voronoi_regions_pixel_based!(p)
    global game_points, canvas_size
    
    # Falls keine Spielpunkte vorhanden sind, Funktion beenden
    if isempty(game_points)
        return
    end

    resolution = 250
    step = canvas_size / resolution # Berechne Schrittweite basierend auf Auflösung
    
    x_coords = 0:step:canvas_size # Erzeuge x-Koordinaten für das Raster
    y_coords = 0:step:canvas_size # Erzeuge y-Koordinaten für das Raster
    z_data = zeros(Int, length(y_coords), length(x_coords)) # Initialisiere Matrix für Spielerbesitzwerte

     # Für jeden Rasterpunkt bestimme den nächstgelegenen Spielpunkt und weise den Besitz zu
    for (i, y) in enumerate(y_coords)
        for (j, x) in enumerate(x_coords)
            test_point = Point(x, y)
            closest_gp = find_closest_game_point(test_point)
            if closest_gp !== nothing
                z_data[i, j] = closest_gp.owner
            end
        end
    end

    custom_colors = cgrad([RGBA(0,0,0,0), RGBA(1,0.5,0.5,0.4), RGBA(0.5,0.5,1,0.4)], [0, 0.5, 1], categorical=true)
 #Zeichne Heatmap auf das Plot-Objekt basierend auf den Besitzdaten
    Plots.heatmap!(p, x_coords, y_coords, z_data,
        color=custom_colors,
        colorbar=false,
        label=""
    )
end

function draw_voronoi_boundaries!(p, voronoi::VoronoiDiagram)
    for (p1, p2) in voronoi.edges # Für jede Kante im Voronoi-Diagramm
        clipped_line = clip_line_to_canvas(p1, p2)#unter：# Schneide die Linie auf die Grenzen der Zeichenfläche zu
        if clipped_line !== nothing
            start_pt, end_pt = clipped_line
            Plots.plot!(p, 
                [start_pt.x, end_pt.x],
                [start_pt.y, end_pt.y],
                color=:black,
                linewidth=1.5,
                alpha=0.8,
                label="")
        end
    end
end

function clip_line_to_canvas(p1::Point, p2::Point)
    # Koordinaten der beiden Punkte extrahieren
    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y
# Grenzen des Zeichenbereichs definieren
    min_x, max_x = 0.0, canvas_size
    min_y, max_y = 0.0, canvas_size
# Zustandscodes für Position relativ zur Zeichenfläche
    INSIDE, LEFT, RIGHT, BOTTOM, TOP = 0, 1, 2, 4, 8

    function compute_outcode(x, y)
        code = INSIDE
        if x < min_x; code |= LEFT; elseif x > max_x; code |= RIGHT; end # Punkt links,rechts außerhalb
        if y < min_y; code |= BOTTOM; elseif y > max_y; code |= TOP; end # Punkt unten,oben außerhalb
        return code
    end
 # Outcodes für beide Endpunkte berechnen
    outcode1 = compute_outcode(x1, y1)
    outcode2 = compute_outcode(x2, y2)
  # Hauptschleife zur Linienbeschneidung
    while true
        # Fall 1: Beide Punkte innerhalb, Linie komplett sichtbar
        if (outcode1 == INSIDE && outcode2 == INSIDE) 
            return (Point(x1, y1), Point(x2, y2))
         # Fall 2: Gemeinsamer Außencode, Linie komplett außerhalb
        elseif (outcode1 & outcode2) != 0
            return nothing
        else
            # Fall 3: Linie muss beschnitten werden
            x_intersect, y_intersect = 0.0, 0.0
            outcode_out = outcode1 != INSIDE ? outcode1 : outcode2

# Berechnung des Schnittpunkts mit der entsprechenden Begrenzung
            if (outcode_out & TOP) != 0
                x_intersect = y1 == y2 ? x1 : x1 + (x2 - x1) * (max_y - y1) / (y2 - y1)
                y_intersect = max_y
            elseif (outcode_out & BOTTOM) != 0
                x_intersect = y1 == y2 ? x1 : x1 + (x2 - x1) * (min_y - y1) / (y2 - y1)
                y_intersect = min_y
            elseif (outcode_out & RIGHT) != 0
                y_intersect = x1 == x2 ? y1 : y1 + (y2 - y1) * (max_x - x1) / (x2 - x1)
                x_intersect = max_x
            elseif (outcode_out & LEFT) != 0
                y_intersect = x1 == x2 ? y1 : y1 + (y2 - y1) * (min_x - x1) / (x2 - x1)
                x_intersect = min_x
            end

            if outcode_out == outcode1
                x1, y1 = x_intersect, y_intersect
                outcode1 = compute_outcode(x1, y1)
            else
                x2, y2 = x_intersect, y_intersect
                outcode2 = compute_outcode(x2, y2)
            end
        end
    end
end

function draw_player_points!(p)
    global game_points
    
    player1_points = [gp for gp in game_points if gp.owner == 1]
    player2_points = [gp for gp in game_points if gp.owner == 2]
    
    if !isempty(player1_points)
        x1 = [gp.point.x for gp in player1_points]
        y1 = [gp.point.y for gp in player1_points]
        Plots.scatter!(p, x1, y1, 
            color=:red, 
            markersize=8, 
            markerstrokewidth=2,
            markerstrokecolor=:darkred,
            label="Player 1 (Rot)")
    end
    
    if !isempty(player2_points)
        x2 = [gp.point.x for gp in player2_points]
        y2 = [gp.point.y for gp in player2_points]
        Plots.scatter!(p, x2, y2, 
            color=:blue, 
            markersize=8, 
            markerstrokewidth=2,
            markerstrokecolor=:darkblue,
            label="Player 2 (Blau)")
    end
end

function get_game_title()
    global current_player, turn_count, game_over, winner
    
    area1, area2 = calculate_areas_voronoi()
    total_area = area1 + area2
    
    title_text = if game_over
        if winner === nothing
            "Spiel beendet - Unentschieden!"
        else
            "Spiel beendet - Player $winner siegt."
        end
    else
        "Zug von Spieler $(current_player)  | Runde $(turn_count + 1)/$max_turns "
    end
    
    if total_area > 0
        pct1 = round(area1 / total_area * 100, digits=1)
        pct2 = round(area2 / total_area * 100, digits=1)
        title_text *= "\nFläche: Rot $(pct1)% | Blau $(pct2)%"
    end
    
    return title_text
end

# ================Interaktion==========================

function click(x::Real, y::Real)
    add_point(x, y)
end

function show()
    show_game()
end

function restart()
    reset_game()
    try
        show_game()
    catch e
        println("false ", e)
    end
end

# ==========================================

println("=== Voronoi-Territoriumsspiel (Wiederhergestellte Version) ===")
println("\nAnleitung:")
println("- `click(x, y)`: Setze einen Punkt an Position (x, y) (Bereich 0–100)")
println("- `show()`: Zeige den aktuellen Spielstand")
println("- `restart()`: Starte das Spiel neu")
println("- Klicke in die Nähe eines vorhandenen Punkts, um ihn zu löschen")

reset_game()

println("\nSpiel gestartet! Bitte beginne. Beispiel: `click(25, 30)`")

